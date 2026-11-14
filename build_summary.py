# FireHistory Atlas: build summary.json from NASA FIRMS MODIS + VIIRS CSV files
# Put the MODIS and VIIRS CSV files in the same folder as this script (in Colab: the "Files" panel).
import glob, json
import numpy as np
import pandas as pd

CELL = 0.25  # grid cell size in degrees (about 25 km)

# ---------- 1. Load and standardize ----------
frames = []
for f in sorted(glob.glob("*.csv")):
    try:
        d = pd.read_csv(f, dtype={"acq_time": str})
    except Exception as e:
        print("Could not read", f, "->", e); continue
    print("File:", f, "| rows:", len(d), "| columns:", list(d.columns))
    if "acq_date" not in d.columns or "latitude" not in d.columns:
        print("  -> skipped (not a FIRMS fire file)"); continue
    inst = str(d["instrument"].iloc[0]).upper() if "instrument" in d.columns and len(d) else ""
    is_viirs = ("VIIRS" in inst) or ("bright_ti4" in d.columns) or ("SV-C2" in f) or ("VIIRS" in f.upper())
    if is_viirs:
        sensor = "VIIRS"
        if "satellite" in d.columns:
            print("  satellites in file:", d["satellite"].astype(str).value_counts().to_dict())
            sel = d[d["satellite"].astype(str).astype(str).str.upper().isin(["N", "NPP", "SNPP", "SUOMI NPP", "S-NPP"])]   # Suomi NPP
            d = sel if len(sel) else d
        cn = pd.to_numeric(d["confidence"], errors="coerce")
        if cn.notna().mean() > 0.5:      # numeric confidence
            d["conf_ok"] = cn >= 30
        else:                            # letters: l / n / h
            d["conf_ok"] = d["confidence"].astype(str).str.lower().isin(["n", "h", "nominal", "high"])
    else:
        sensor = "MODIS"
        d["conf_ok"] = pd.to_numeric(d["confidence"], errors="coerce") >= 30
    d = d.dropna(subset=["latitude", "longitude", "acq_date", "acq_time", "frp"])
    print("  -> treated as", sensor, "| dates:", d["acq_date"].min(), "to", d["acq_date"].max(), "| usable rows:", len(d))
    d["sensor"] = sensor
    d["acq_time"] = d["acq_time"].astype(str).str.zfill(4)
    d["dt_utc"] = pd.to_datetime(d["acq_date"] + " " + d["acq_time"], format="%Y-%m-%d %H%M", utc=True)
    cols = ["sensor", "latitude", "longitude", "acq_date", "dt_utc", "frp", "conf_ok"]
    if "type" in d.columns:
        cols.append("type")
    frames.append(d[cols])
raw = pd.concat(frames, ignore_index=True)
if "type" not in raw.columns:
    raw["type"] = 0
sensors_found = sorted(raw["sensor"].unique())
print("Sensors found:", sensors_found)
if sensors_found != ["MODIS", "VIIRS"]:
    raise SystemExit("Need BOTH a MODIS file and a VIIRS file in the folder. Found: %s" % sensors_found)

# Check for cut-off uploads
for sname in sensors_found:
    dd = raw[raw["sensor"] == sname]["acq_date"]
    span = (pd.to_datetime(dd.max()) - pd.to_datetime(dd.min())).days
    if span < 365 * 3:
        raise SystemExit("STOP: the %s file only covers %s to %s (%d rows). The upload looks incomplete. "
                         "Delete it in the file panel and upload it again, waiting until it finishes."
                         % (sname, dd.min(), dd.max(), len(dd)))

# Raw counts per year BEFORE filtering: this shows the sensor mismatch
raw["year"] = raw["acq_date"].str[:4]
raw_by_year = raw.groupby(["sensor", "year"]).size().unstack(0).fillna(0).astype(int)
print("\nRaw detections per year (note how different the two sensors are):")
print(raw_by_year)

# Keep vegetation-type fires with nominal/high confidence
d = raw[(raw["type"] == 0) & raw["conf_ok"]].copy()

# Only use dates where both sensors have data
start = max(d[d.sensor == s]["acq_date"].min() for s in sensors_found)
end = min(d[d.sensor == s]["acq_date"].max() for s in sensors_found)
d = d[(d["acq_date"] >= start) & (d["acq_date"] <= end)]

# ---------- 2. Grid cells ----------
d["cx"] = np.floor(d["longitude"] / CELL).astype(int)
d["cy"] = np.floor(d["latitude"] / CELL).astype(int)
d["cell"] = (d["cy"] + 2000) * 10000 + (d["cx"] + 2000)

# ---------- 3. Daily metrics per sensor ----------
days = pd.date_range(start, end, freq="D")
daily = d.groupby(["sensor", "acq_date"]).agg(detections=("frp", "size"), frp=("frp", "sum"),
                                              cells=("cell", "nunique")).reset_index()
daily["acq_date"] = pd.to_datetime(daily["acq_date"])
def series(sensor, col):
    s = daily[daily.sensor == sensor].set_index("acq_date")[col]
    return s.reindex(days, fill_value=0).astype(float)
M = {c: series("MODIS", c) for c in ["cells", "frp", "detections"]}
V = {c: series("VIIRS", c) for c in ["cells", "frp", "detections"]}

# ---------- 4. Overlap ratio (VIIRS / MODIS) by month ----------
mon = pd.DataFrame({"m_cells": M["cells"], "v_cells": V["cells"], "m_frp": M["frp"], "v_frp": V["frp"],
                    "m_det": M["detections"], "v_det": V["detections"]}).groupby(lambda x: x.strftime("%Y-%m")).sum()
mon["cells_ratio"] = mon["v_cells"] / mon["m_cells"].replace(0, np.nan)
mon["frp_ratio"] = mon["v_frp"] / mon["m_frp"].replace(0, np.nan)
mon["count_ratio"] = mon["v_det"] / mon["m_det"].replace(0, np.nan)
def summ(x):
    x = x.dropna()
    return {"median": round(float(x.median()), 2), "p25": round(float(x.quantile(.25)), 2),
            "p75": round(float(x.quantile(.75)), 2)}
ratio_summary = {"cell_days": summ(mon["cells_ratio"]), "frp": summ(mon["frp_ratio"]),
                 "raw_counts": summ(mon["count_ratio"])}
print("\nVIIRS/MODIS ratio (monthly median and spread):")
print(json.dumps(ratio_summary, indent=1))

# ---------- 5. Unusual days (median and MAD, same time of year in OTHER years) ----------
def robust_z(s):
    doy = s.index.dayofyear.values.clip(max=365); yr = s.index.year.values; v = s.values
    out = np.zeros(len(v))
    for i in range(len(v)):
        dist = np.abs(doy - doy[i]); dist = np.minimum(dist, 365 - dist)
        other = v[(yr != yr[i]) & (dist <= 7)]
        if len(other) < 5:
            continue
        med = np.median(other); mad = np.median(np.abs(other - med))
        out[i] = (v[i] - med) / max(1.4826 * mad, 1.0)
    return out
zm, zv = robust_z(M["cells"]), robust_z(V["cells"])
zc = (zm + zv) / 2
THRESH = 3.0

flag = np.where(zc >= THRESH)[0]
episodes, cur = [], []
for i in flag:
    if cur and i - cur[-1] > 2:
        episodes.append(cur); cur = []
    cur.append(i)
if cur:
    episodes.append(cur)
ep_list = []
for e in episodes:
    p = e[int(np.argmax(zc[e]))]
    ep_list.append({"start": days[e[0]].strftime("%Y-%m-%d"), "end": days[e[-1]].strftime("%Y-%m-%d"),
                    "peak_date": days[p].strftime("%Y-%m-%d"), "peak_score": round(float(zc[p]), 1),
                    "flagged_days": len(e), "total_score": round(float(zc[e].sum()), 1), "both_sensors_agree": bool(zm[p] >= THRESH and zv[p] >= THRESH),
                    "modis_cells": int(M["cells"].iloc[p]), "viirs_cells": int(V["cells"].iloc[p])})
ep_list = sorted(ep_list, key=lambda x: -x["total_score"])[:20]   # biggest = strongest AND longest
print("\nTop unusual periods:")
for e in ep_list[:10]:
    print(e["start"], "to", e["end"], "| peak", e["peak_date"], "score", e["peak_score"], "| days", e["flagged_days"], "| sensors agree" if e["both_sensors_agree"] else "| one sensor")

# ---------- 6. Map data ----------
d["month"] = d["acq_date"].str[:7]
cm = d.groupby(["month", "cell", "cx", "cy", "sensor"])["acq_date"].nunique().unstack("sensor").fillna(0).reset_index()
cell_month = {}
for r in cm.itertuples(index=False):
    cell_month.setdefault(r.month, []).append([round((r.cy + .5) * CELL, 3), round((r.cx + .5) * CELL, 3),
                                              int(getattr(r, "MODIS", 0)), int(getattr(r, "VIIRS", 0))])
anom_cells = {}
for e in ep_list:
    sub = d[d["acq_date"] == e["peak_date"]]
    g = sub.groupby(["cx", "cy"]).agg(m=("sensor", lambda s: int((s == "MODIS").any())),
                                      v=("sensor", lambda s: int((s == "VIIRS").any())),
                                      frp=("frp", "sum")).reset_index()
    anom_cells[e["peak_date"]] = [[round((r.cy + .5) * CELL, 3), round((r.cx + .5) * CELL, 3), r.m, r.v, round(float(r.frp), 1)]
                                  for r in g.itertuples(index=False)]

# ---------- 7. Write summary.json ----------
r = lambda a: [round(float(x), 1) for x in a]
out = {
    "meta": {
        "region": "California area (bounding box %.1f to %.1f N, %.1f to %.1f W)" % (
            raw["latitude"].min(), raw["latitude"].max(), -raw["longitude"].max(), -raw["longitude"].min()),
        "start": start, "end": end, "grid_degrees": CELL, "anomaly_threshold": THRESH,
        "sources": "NASA FIRMS MODIS Collection 6.1 (Terra and Aqua) and VIIRS S-NPP, Standard Processing",
        "filters": "Vegetation-fire type only; MODIS confidence 30 or higher; VIIRS nominal or high confidence",
        "notice": "Hotspots are satellite detections of thermal anomalies, not confirmed fires and not predictions."},
    "daily": {"start": start, "modis_cells": [int(x) for x in M["cells"]], "viirs_cells": [int(x) for x in V["cells"]],
              "modis_frp": r(M["frp"]), "viirs_frp": r(V["frp"]), "z_modis": r(zm), "z_viirs": r(zv), "z_combined": r(zc)},
    "monthly_ratio": [{"month": k, "cells": None if pd.isna(v.cells_ratio) else round(v.cells_ratio, 2),
                       "frp": None if pd.isna(v.frp_ratio) else round(v.frp_ratio, 2),
                       "counts": None if pd.isna(v.count_ratio) else round(v.count_ratio, 2)}
                      for k, v in mon.iterrows()],
    "ratio_summary": ratio_summary,
    "raw_by_year": {s: {y: int(raw_by_year.loc[y, s]) for y in raw_by_year.index} for s in sensors_found},
    "episodes": ep_list, "cell_month": cell_month, "anomaly_cells": anom_cells}
with open("summary.json", "w") as f:
    json.dump(out, f, separators=(",", ":"))
import os
print("\nDone. summary.json size: %.0f KB" % (os.path.getsize("summary.json") / 1024))
