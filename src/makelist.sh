
#!/bin/bash

find . -type f | grep .cpp | xargs -I{} basename {} '.cpp' | paste -sd ' ' > libcapdddes.txt

find . -type d | paste -sd ' ' > libvpath.txt