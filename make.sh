#!/bin/bash

g++ main.cc \
    -I/Users/changhyonha/Software/cry_v1.7/src \
    -L/Users/changhyonha/Software/cry_v1.7/lib \
    -lCRY \
    $(root-config --cflags --libs) \
    -o hawlu_inmotion
