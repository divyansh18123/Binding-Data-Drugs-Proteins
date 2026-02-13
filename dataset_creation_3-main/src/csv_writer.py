import threading
import pandas as pd

_csv_lock = threading.Lock()

def safe_to_csv(df: pd.DataFrame, path: str, header: bool, mode: str):
    with _csv_lock:
        df.to_csv(path, index=False, header=header, mode=mode)
