#!/bin/sh
# Rebuild motif_engine.{js,wasm} from src/. Requires Emscripten (emcc).
set -eu
cd "$(dirname "$0")"
emcc -std=c11 -O2 -Wall -Wextra -D_DEFAULT_SOURCE src/*.c \
  -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createMotifEngine \
  -sENVIRONMENT=web,node -sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=4GB \
  -sEXPORTED_FUNCTIONS=_motif_compare,_motif_malloc,_motif_free \
  -sEXPORTED_RUNTIME_METHODS=UTF8ToString,stringToNewUTF8,HEAPF32 \
  -o motif_engine.js
ls -l motif_engine.js motif_engine.wasm
