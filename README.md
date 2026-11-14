# FireHistory Atlas

**A harmonized burning-activity calendar built from NASA MODIS and VIIRS active-fire data.**

NASA Space Apps Challenge 2026 · Challenge: Harmonization of MODIS and VIIRS Hot Spots

🔗 **Live demo:** https://anaha06.github.io/firehistory-atlas/

## What it does

FireHistory Atlas lets users explore past fire activity for [your region] from [start year] to [end year]. It shows a calendar of daily activity, a map of detections, a ranked list of unusually active days, and a plain-language explanation of how the two satellite sensors differ.

## Why it matters

MODIS and VIIRS detect fires differently (VIIRS has finer pixels and detects smaller fires), so their raw counts cannot be compared directly. This project uses comparable measures so users can study historical patterns more fairly.

## How we handle the sensor differences

- We do not convert individual pixels between sensors.
- We count **unique grid-cell days** with at least one detection instead of raw detections.
- We also sum **Fire Radiative Power (FRP)** per day.
- We compare the two sensors over the period where both operate, and report the VIIRS/MODIS ratio.
- Unusual days are flagged against the typical level for the same time of year.

## What this is not

Hotspots are satellite detections of thermal anomalies, not confirmed fires and not predictions. A day with no detections does not prove there was no fire (clouds and satellite overpass times matter). This tool does not predict future fires and does not replace official emergency alerts.

## Data

- NASA FIRMS (Fire Information for Resource Management System): MODIS Collection 6.1 and VIIRS S-NPP, Standard Processing archive

## Built with

Python, Pandas, Leaflet, d3, HTML/JavaScript

## Use of AI

AI tools were used to help plan the project and write parts of the code. [Describe what you changed or wrote yourself.]

