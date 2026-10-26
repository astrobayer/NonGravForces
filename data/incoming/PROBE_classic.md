# fetchData probe report, 2026-10-26 18:28 UTC

Written by `tools/fetchData.py` on a CI runner (hosted: GitHub; classic: self-hosted). Staging only: see the script's docstring. Each line is a request and what came back.

## Classic GRACE Level-1B, 2008-01-01 (N-08)

- `https://isdc-data.gfz.de/grace/Level-1B/` -> failed -- URLError: <urlopen error [WinError 10060] Ein Verbindungsversuch ist fehlgeschlagen, da die Gegenstelle nach einer bestimmten Zeitspanne nicht richtig reagiert hat, oder die hergestellte Verbindung war fehlerhaft, da der verbundene Host nicht reagiert hat>
- `ftp://isdcftp.gfz-potsdam.de/grace/Level-1B/` -> failed -- URLError: <urlopen error [Errno 11001] getaddrinfo failed>
`ftp://isdcftp.gfz-potsdam.de/grace/Level-1B/`: 
