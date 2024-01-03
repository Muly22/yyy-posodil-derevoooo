main: main.c
	mkdir -p save
	touch save/save.txt
	gcc -o main main.c -lm
	x86_64-w64-mingw32-gcc -o main.exe main.c -lm
clean:
	rm -fr main main.exe save/save.txt