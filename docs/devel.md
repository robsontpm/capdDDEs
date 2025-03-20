Dev's notes
===========

Notes for developers, everything we need to remeber for future, and
whta is done very rarely. 


Submodules
----------

Currently the repository has those submodules:
 
 - CAPD

In the case you want to list submodules, just write:

```git submodule```

To add the CAPD submodule write

```git submodule add https://github.com/CAPDGroup/CAPD.git external/capd```

to update submodule in your cloned repostiroy write:

If this is freshly copied repository: 

```
git submodule init
git submodule update
```
Then, everytime you want to get the update from submodule:

```git submodule update --remote```

This should work most of the time. 






