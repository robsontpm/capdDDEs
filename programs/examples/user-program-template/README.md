A template for a user program
=============================

Please copy this whole folder to ```programs/devel/myprogramname``` (put your favourity
name in place of ```myprogramname```. You should have:

```
programs/devel/myprogramname
 - Makefile
 - program.cpp
 - README.md
 - setup.h
 - setup.cpp
```

You might alse rename other files to what you like, but you need to modify the ```Makefile``` accordingly. 

Then, inside the ```programs/devel/myprogramname```
you can invoke compilation by 

```
make all
```

to cleanup

```
make clean
```

The executables will go to ```programs/devel/myprogramname/bin/```.

```
cd bin
./program
```

Now, modify ```setup.h``` for your equation and 
then ```program.cpp``` to do your computations!
    