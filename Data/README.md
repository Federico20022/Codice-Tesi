# Data

This folder is empty on purpose. The SPX option quotes come from LSEG Workspace (Refinitiv) and are under licence, so they are not distributed.

Put here one file per date, named `AAAAMMGGSPX_merged.csv`, with the columns

```
Instrument,Root,PutCall,CF_LAST,CF_BID,CF_ASK,STRIKE_PRC,EXPIR_DATE
```

and one row with `Root = SPOT`. See the main README for the run order.
