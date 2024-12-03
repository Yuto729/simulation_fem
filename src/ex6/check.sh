#!/bin/bash

set -e
gcc check.c -o check -lm

./check
