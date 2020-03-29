# WSL error
```
: not foundele_profiler.sh: 2: ./DAJIN/allele_profiler.sh:
: not foundele_profiler.sh: 6: ./DAJIN/allele_profiler.sh:
./DAJIN/allele_profiler.sh: 7: set: Illegal option -
```
This error is caused by `CRLF`, a bytecode that can be used to mark a line break in a text file. Please open allele_profiler.sh by your text editor (e.g. Visual Studio Code) and change `CRLF` to `LF`.
