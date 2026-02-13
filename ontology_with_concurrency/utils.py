# src/http_aio_wrapper_fast.py
import asyncio
import atexit
import random
import threading
import logging
from typing import Any, Dict, Optional, Callable, Iterable
import aiohttp

# optional faster json parser
try:
    import orjson as _json
    _loads = _json.loads
except Exception:
    import json as _stdjson
    _loads = lambda b: _stdjson.loads(b.decode("utf-8"))

# optional uvloop for faster event loop (unix only)
try:
    import uvloop
    _use_uvloop = True
except Exception:
    _use_uvloop = False

# optional small LRU cache
try:
    from cachetools import LRUCache, cached
    _have_cache = True
except Exception:
    _have_cache = False

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

DEFAULT_HEADERS = {
    "Accept": "application/json, text/plain, */*",
    "Accept-Encoding": "gzip, deflate, br",
    "User-Agent": "MyApp/AsyncHTTP/fast/1.0",
}

class _AioServiceFast:
    def __init__(
        self,
        headers: Optional[Dict[str, str]] = None,
        concurrency: int = 100,
        timeout: float = 10.0,
        retries: int = 3,
        backoff_factor: float = 0.25,
        enable_cache: bool = False,
        cache_maxsize: int = 1024,
    ):
        """
        Faster aiohttp wrapper with orjson + uvloop + optional LRU cache.

        - concurrency: semaphore limit (how many in-flight requests you'll allow)
        - timeout: per-request total timeout
        - retries/backoff_factor: for 429/5xx or network errors
        - enable_cache: if True, caches JSON responses by URL (simple LRU)
        """
        self._headers = headers or DEFAULT_HEADERS.copy()
        self._concurrency = concurrency
        self._timeout = timeout
        self._retries = max(0, retries)
        self._backoff_factor = backoff_factor
        self._enable_cache = enable_cache and _have_cache

        # background loop/session
        self._loop: Optional[asyncio.AbstractEventLoop] = None
        self._thread: Optional[threading.Thread] = None
        self._session: Optional[aiohttp.ClientSession] = None
        self._sem: Optional[asyncio.Semaphore] = None
        self._started = threading.Event()
        self._closed = False

        # small LRU cache for repeated calls (URL key)
        if self._enable_cache:
            self._cache = LRUCache(maxsize=cache_maxsize)
        else:
            self._cache = None

        self._start_background_loop()

    def _start_background_loop(self):
        if self._thread and self._thread.is_alive():
            return

        def _run():
            # optionally set uvloop (UNIX)
            if _use_uvloop:
                uvloop.install()
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            self._loop = loop

            async def _init():
                # Use a TCPConnector tuned for many connections.
                # We set limit=0 to remove connector-level limit and use the semaphore to control concurrency.
                connector = aiohttp.TCPConnector(limit=0, enable_cleanup_closed=True, ttl_dns_cache=300)
                timeout = aiohttp.ClientTimeout(total=self._timeout)
                # pass trust_env=False to avoid proxy overhead unless you need proxies
                self._session = aiohttp.ClientSession(headers=self._headers, timeout=timeout, connector=connector, trust_env=False)
                self._sem = asyncio.Semaphore(self._concurrency)
                self._started.set()

            loop.run_until_complete(_init())
            try:
                loop.run_forever()
            finally:
                try:
                    loop.run_until_complete(self._shutdown())
                except Exception:
                    pass
                loop.close()

        self._thread = threading.Thread(target=_run, daemon=True)
        self._thread.start()
        if not self._started.wait(timeout=10.0):
            raise RuntimeError("Failed to start background aiohttp service")

    async def _shutdown(self):
        if self._session and not self._session.closed:
            try:
                await self._session.close()
            except Exception:
                pass
        self._session = None

    def close(self):
        if self._closed:
            return
        self._closed = True
        if not self._loop:
            return
        # stop loop safely
        def _stop():
            if self._loop and self._loop.is_running():
                self._loop.stop()
        self._loop.call_soon_threadsafe(_stop)
        if self._thread:
            self._thread.join(timeout=5.0)

    def _schedule(self, coro: asyncio.coroutines) -> "asyncio.Future":
        if self._closed:
            raise RuntimeError("Service already closed")
        if not self._loop:
            raise RuntimeError("Background loop not available")
        return asyncio.run_coroutine_threadsafe(coro, self._loop)

    async def _fetch_coro(
        self,
        url: str,
        *,
        params: Optional[Dict[str, str]] = None,
        headers: Optional[Dict[str, str]] = None,
        expect_json: bool = True,
    ) -> Optional[Any]:
        """core coroutine: uses resp.read() + orjson loads for speed"""
        if self._session is None:
            raise RuntimeError("Session not ready")
        session = self._session
        sem = self._sem  # type: ignore

        # merge headers: use shallow copy of session default headers fallback
        req_headers = dict(session._default_headers) if hasattr(session, "_default_headers") else dict(self._headers)
        if headers:
            req_headers.update(headers)

        # simple URL-only cache
        cache_key = url if self._enable_cache else None
        if self._enable_cache and cache_key in self._cache:
            return self._cache[cache_key]

        attempt = 0
        while True:
            try:
                async with sem:
                    async with session.get(url, params=params, headers=req_headers) as resp:
                        status = resp.status
                        if 200 <= status < 300:
                            # read bytes then parse with orjson (fast)
                            b = await resp.read()
                            if expect_json:
                                try:
                                    obj = _loads(b)
                                except Exception:
                                    # fallback to text
                                    obj = b.decode("utf-8", errors="replace")
                                if cache_key is not None:
                                    self._cache[cache_key] = obj
                                return obj
                            else:
                                txt = b.decode("utf-8", errors="replace")
                                if cache_key is not None:
                                    self._cache[cache_key] = txt
                                return txt

                        # retry on 429 or 5xx
                        if status == 429 or 500 <= status < 600:
                            if attempt < self._retries:
                                backoff = self._backoff_factor * (2 ** attempt)
                                jitter = random.uniform(0, 0.1 * backoff)
                                await asyncio.sleep(backoff + jitter)
                                attempt += 1
                                continue
                            else:
                                logger.debug("Http status %s for %s", status, url)
                                return None

                        # other client errors: return None
                        return None

            except (aiohttp.ClientError, asyncio.TimeoutError) as exc:
                if attempt < self._retries:
                    backoff = self._backoff_factor * (2 ** attempt)
                    jitter = random.uniform(0, 0.1 * backoff)
                    await asyncio.sleep(backoff + jitter)
                    attempt += 1
                    continue
                logger.debug("Network error fetching %s: %s", url, exc)
                return None

    # sync wrapper to call from normal synchronous code
    def fetch(self, url: str, params: Optional[Dict[str, str]] = None, headers: Optional[Dict[str, str]] = None, expect_json: bool = True, timeout: Optional[float] = None) -> Optional[Any]:
        # ensure not called from inside a running loop in this thread
        try:
            _ = asyncio.get_running_loop()
            raise RuntimeError("fetch() must not be called inside an active asyncio event loop; use fetch_async_* instead")
        except RuntimeError:
            pass

        fut = self._schedule(self._fetch_coro(url, params=params, headers=headers, expect_json=expect_json))
        wait_timeout = (timeout or self._timeout) + 2.0
        try:
            return fut.result(timeout=wait_timeout)
        except Exception as e:
            logger.debug("Exception waiting for fetch result %s: %s", url, e)
            try:
                fut.cancel()
            except Exception:
                pass
            return None

    # async interface for code that is already async
    async def fetch_async(self, url: str, params: Optional[Dict[str, str]] = None, headers: Optional[Dict[str, str]] = None, expect_json: bool = True) -> Optional[Any]:
        # If called from inside the same loop used by the service, call directly
        if asyncio.get_event_loop() is self._loop:
            return await self._fetch_coro(url, params=params, headers=headers, expect_json=expect_json)
        fut = self._schedule(self._fetch_coro(url, params=params, headers=headers, expect_json=expect_json))
        return fut.result()

    async def fetch_many_async(self, urls: Iterable[str], *, batch_size: int = 200) -> Dict[str, Optional[Any]]:
        """
        Async batch fetch using gather in chunks - efficient and low-overhead.
        """
        out = {}
        url_list = list(urls)
        for i in range(0, len(url_list), batch_size):
            chunk = url_list[i:i+batch_size]
            tasks = [asyncio.create_task(self._fetch_coro(u, expect_json=True)) for u in chunk]
            # gather ensures tasks run concurrently and we avoid as_completed overhead
            results = await asyncio.gather(*tasks, return_exceptions=True)
            for u, r in zip(chunk, results):
                out[u] = None if isinstance(r, Exception) else r
        return out

# singleton instance
_service = _AioServiceFast(concurrency=80, timeout=8.0, retries=2, backoff_factor=0.25, enable_cache=False)

atexit.register(lambda: _service.close())

# public API (same names as before)
def fetch_json(url: str, headers: Optional[Dict[str, str]] = None, params: Optional[Dict[str, str]] = None, timeout: Optional[float] = None) -> Optional[Any]:
    return _service.fetch(url, params=params, headers=headers, expect_json=True, timeout=timeout)

def fetch_text(url: str, headers: Optional[Dict[str, str]] = None, params: Optional[Dict[str, str]] = None, timeout: Optional[float] = None) -> Optional[str]:
    return _service.fetch(url, params=params, headers=headers, expect_json=False, timeout=timeout)

async def fetch_json_async(url: str, headers: Optional[Dict[str, str]] = None, params: Optional[Dict[str, str]] = None) -> Optional[Any]:
    return await _service.fetch_async(url, params=params, headers=headers, expect_json=True)

async def fetch_text_async(url: str, headers: Optional[Dict[str, str]] = None, params: Optional[Dict[str, str]] = None) -> Optional[str]:
    return await _service.fetch_async(url, params=params, headers=headers, expect_json=False)

async def fetch_many_async(urls: Iterable[str], *, batch_size: int = 200) -> Dict[str, Optional[Any]]:
    return await _service.fetch_many_async(urls, batch_size=batch_size)
