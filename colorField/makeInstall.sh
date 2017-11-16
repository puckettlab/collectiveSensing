#!/bin/bash
#for my program
#g++ my_program.cc -o my_program `freetype-config --cflags`  -I/media/data1/aawaz/create_field/include  -L/media/data1/aawaz/create_field/lib -lpngwriter -lpng  -lz -lfreetype


#for create_field
g++ main.cpp -o my_program `freetype-config --cflags` -static -I/media/data1/aawaz/code/cpp/create_field/include  -L/media/data1/aawaz/code/cpp/create_field/lib -I/usr/include/freetype2 -lpngwriter -lpng  -lz -lfreetype -lgsl
