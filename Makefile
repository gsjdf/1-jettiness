#!/bin/sh




all: fig1a fig1
new: fig1an fig1n
full:fig1afull fig1full
fig1a:
	create-nlojet-user fig1a.cc alew.f -Wx,-O3,-Wall,-I/usr/local/include/LHAPDF ,-I/usr/local/include/fastjet -Wf,-O3,-Wall -Wl,-lfastjet,-lLHAPDF -o fig1a
	
fig1an:
	create-nlojet-user fig1a_new.cc alew.f -Wx,-O3,-Wall,-I/usr/local/include/LHAPDF ,-I/usr/local/include/fastjet -Wf,-O3,-Wall -Wl,-lfastjet,-lLHAPDF -o fig1a_new

fig1afull:
	create-nlojet-user fig1a_new_full.cc alew.f -Wx,-O3,-Wall,-I/usr/local/include/LHAPDF ,-I/usr/local/include/fastjet -Wf,-O3,-Wall -Wl,-lfastjet,-lLHAPDF -o fig1a_full
	
fig1:
	create-nlojet-user fig1.cc alew.f -Wx,-O3,-Wall,-I/usr/local/include/LHAPDF ,-I/usr/local/include/fastjet -Wf,-O3,-Wall -Wl,-lfastjet,-lLHAPDF -o fig1
fig1n:
	create-nlojet-user fig1_new.cc alew.f -Wx,-O3,-Wall,-I/usr/local/include/LHAPDF ,-I/usr/local/include/fastjet -Wf,-O3,-Wall -Wl,-lfastjet,-lLHAPDF -o fig1_new
fig1full:
	create-nlojet-user fig1_new_full.cc alew.f -Wx,-O3,-Wall,-I/usr/local/include/LHAPDF ,-I/usr/local/include/fastjet -Wf,-O3,-Wall -Wl,-lfastjet,-lLHAPDF -o fig1_new_full

clean:
	rm -rf .libs .obj *.la
	rm -rf *.o *~
