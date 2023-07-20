#!/bin/bash

tail -n+3 *score.sc | sort -k2 -n | head -n 1
