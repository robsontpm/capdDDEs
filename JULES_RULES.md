- When creating a new feature branch for working on a task, make sure it's name starts with 'jules', e.g. 'jules/my-new-feature'
- Never attempt cmake . in the root directory of the program!
- Always build, make and compile from within ./build/ directory!
- Do not remove anything from .gitignore when not asked to!
- Do not alter './external/capd-build.sh', if you need to change something, use a copy of the 'capd-build.sh', e.g. 'jules-capd-build.sh'
- IMPORTANT: For test coverage use lconv version 1.16, as the never versions has problems with templates in C++.





