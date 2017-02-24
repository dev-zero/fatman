#!/bin/sh

# dropping the local database first, therefore have to connect to different db
ssh \
    tctdb /usr/pgsql-9.5/bin/pg_dump --clean --create fatman | psql postgres

# copy the file data
rsync -av tctdb:/var/lib/fatman/data/results/ /tmp/resultdata/
