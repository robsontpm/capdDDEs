GROUND RULES
------------

Those rules should be followed to the letter:

- When creating a new feature branch for working on a task, make sure it's name starts with `jules`, e.g. `jules/my-new-feature`
- Never attempt `cmake .` in the root directory of the program!
- Always build, make and compile from within `./build/ directory`!
- Do not remove anything from `.gitignore` when not asked to!
- Do not alter `./external/capd-build.sh`, if you need to change something, use a copy of the `capd-build.sh`, e.g. `jules-capd-build.sh`
- If you need some kind of dependency, install them in `${HOME}/deps` and add `${HOME}/deps` to `${PATH}` by `export PATH=$PATH:/path/to/lcov` in your environment. Do not alter the `.sh` files with specific locations / exports / etc. Keep the build files clean. 
- Use `lcov` version 1.16 for testing code coverage! 





