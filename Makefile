all:
	odin build . -out:fastDTW.exe --debug

opti:
	# odin build . -out:fastDTW.exe -o:speed -no-bounds-check -disable-assert
	odin build . -out:fastDTW.exe -o:aggressive -microarch:native -no-bounds-check -disable-assert

clean:
	rm fastDTW.exe

run:
	./fastDTW.exe


