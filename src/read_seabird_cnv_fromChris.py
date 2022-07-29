import numpy as np
import pandas 
import bz2
import gzip
import linecache
import re
import warnings
import zipfile
from datetime import datetime
from io import StringIO
from pathlib import Path


def _parse_seabird(lines, ftype="cnv"):
    # Initialize variables.
    lon = lat = time = None, None, None
    skiprows = 0

    metadata = {}
    header, config, names = [], [], []
    for k, line in enumerate(lines):
        line = line.strip()

        # Only cnv has columns names, for bottle files we will use the variable row.
        if ftype == "cnv":
            if "# name" in line:
                name, unit = line.split("=")[1].split(":")
                name, unit = list(map(_normalize_names, (name, unit)))
                names.append(name)

        # Seabird headers starts with *.
        if line.startswith("*"):
            header.append(line)

        # Seabird configuration starts with #.
        if line.startswith("#"):
            config.append(line)

        # NMEA position and time.
        if "NMEA Latitude" in line:
            hemisphere = line[-1]
            lat = line.strip(hemisphere).split("=")[1].strip()
            lat = np.float_(lat.split())
            if hemisphere == "S":
                lat = -(lat[0] + lat[1] / 60.0)
            elif hemisphere == "N":
                lat = lat[0] + lat[1] / 60.0
            else:
                raise ValueError("Latitude not recognized.")
        if "NMEA Longitude" in line:
            hemisphere = line[-1]
            lon = line.strip(hemisphere).split("=")[1].strip()
            lon = np.float_(lon.split())
            if hemisphere == "W":
                lon = -(lon[0] + lon[1] / 60.0)
            elif hemisphere == "E":
                lon = lon[0] + lon[1] / 60.0
            else:
                raise ValueError("Latitude not recognized.")
        if "NMEA UTC (Time)" in line:
            time = line.split("=")[-1].strip()
            # Should use some fuzzy datetime parser to make this more robust.
            time = datetime.strptime(time, "%b %d %Y %H:%M:%S")

        # cnv file header ends with *END* while
        if ftype == "cnv":
            if line == "*END*":
                skiprows = k + 1
                break
        else:  # btl.
            # There is no *END* like in a .cnv file, skip two after header info.
            if not (line.startswith("*") | line.startswith("#")):
                # Fix commonly occurring problem when Sbeox.* exists in the file
                # the name is concatenated to previous parameter
                # example:
                #   CStarAt0Sbeox0Mm/Kg to CStarAt0 Sbeox0Mm/Kg (really two different params)
                line = re.sub(r"(\S)Sbeox", "\\1 Sbeox", line)

                names = line.split()
                skiprows = k + 2
                break
    if ftype == "btl":
        # Capture stat names column.
        names.append("Statistic")
    metadata.update(
        {
            "header": "\n".join(header),
            "config": "\n".join(config),
            "names": names,
            "skiprows": skiprows,
            "time": time,
            "lon": lon,
            "lat": lat,
        }
    )
    return metadata


def from_cnv(fname):
    
    f = _read_file(fname)
    metadata = _parse_seabird(f.readlines(), ftype="cnv")
    
    f.seek(0)
    df = pandas.read_fwf(
        f,
        header=None,
        index_col=None,
        names=metadata["names"],
        skiprows=metadata["skiprows"],
        delim_whitespace=True,
        widths=[11] * len(metadata["names"]),
    )
    f.close()

    key_set = False
    prkeys = ["prDM", "prdM", "pr"]
    for prkey in prkeys:
        try:
            df.set_index(prkey, drop=True, inplace=True)
            key_set = True
        except KeyError:
            continue
    if not key_set:
        raise KeyError(
            f"Could not find pressure field (supported names are {prkeys})."
        )
    df.index.name = "Pressure [dbar]"

    name = _basename(fname)[1]
        

    dtypes = {"bpos": int, "pumps": bool, "flag": bool}
    for column in df.columns:
        if column in dtypes:
            df[column] = df[column].astype(dtypes[column])
#        elif column=='longitude':
#            print(df[column])
#            df[column] = df[column].astype(float)
#            
#            print(df[column])
#            dsa
      
        else:
            try:
                df[column] = df[column].astype(float)
            except ValueError:
                warnings.warn("Could not convert %s to float." % column)

    metadata["name"] = str(name)
    setattr(df, "_metadata", metadata)
    return df


def _read_file(fname):
    if not isinstance(fname, Path):
        fname = Path(fname).resolve()

    extension = fname.suffix.lower()
    if extension in [".gzip", ".gz", ".bz2", ".zip"]:
        contents = _open_compressed(fname)
    elif extension in [".cnv", ".edf", ".txt", ".ros", ".btl"]:
        contents = fname.read_bytes()
    else:
        raise ValueError(
            f"Unrecognized file extension. Expected .cnv, .edf, .txt, .ros, or .btl got {extension}"
        )
    # Read as bytes but we need to return strings for the parsers.
    text = contents.decode(encoding="utf-8", errors="replace")
    return StringIO(text)


def _normalize_names(name):
    name = name.strip()
    name = name.strip("*")
    return name

def _basename(fname):
    """Return file name without path."""
    if not isinstance(fname, Path):
        fname = Path(fname)
        path, name, ext = fname.parent, fname.stem, fname.suffix
        return path, name, ext
      
    else:
        try:
            df[column] = df[column].astype(float)
        except ValueError:
            warnings.warn("Could not convert %s to float." % column)

    metadata["name"] = str(name)
    setattr(df, "_metadata", metadata)
    return df

def _read_file(fname):
    if not isinstance(fname, Path):
        fname = Path(fname).resolve()

    extension = fname.suffix.lower()
    if extension in [".gzip", ".gz", ".bz2", ".zip"]:
        contents = _open_compressed(fname)
    elif extension in [".cnv", ".edf", ".txt", ".ros", ".btl"]:
        contents = fname.read_bytes()
    else:
        raise ValueError(
            f"Unrecognized file extension. Expected .cnv, .edf, .txt, .ros, or .btl got {extension}"
        )
    # Read as bytes but we need to return strings for the parsers.
    text = contents.decode(encoding="utf-8", errors="replace")
    return StringIO(text)


