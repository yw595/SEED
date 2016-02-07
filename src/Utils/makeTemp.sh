#!/bin/bash

cat $1 | tr "&" "_" | grep -o -P "rxn\d+,.*(\d|-)+\.(\d|-)+\.(\d|-)+\.(\d|-)+" > $2
