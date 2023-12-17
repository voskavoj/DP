call "C:\ProgramData\miniconda3\condabin\activate.bat"
activate gnuradioenv
cd "C:\Git\Personal\DP\Data"
iridium-extractor -D 4 "C:\Git\Personal\DP\GNURadio\plutosdr-soapy.conf" > output.bits
pause