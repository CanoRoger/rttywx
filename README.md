# HF RTTY marine weather forecast decoder

This library is a modification of libnavtex by Franco Venturi to decode instead RTTY weather forecast in short wave.
The decoder code is extracted from fldigi source code, reduced and adapted to be used as a separated library.

Later on, this library is used in GNURadio gr-rttyrx.

## How to build and install

```
git clone https://github.com/CanoRoger/librttywx.git
cd librttywx
mkdir build
cd build
cmake ..
make
sudo make install
```

## How to run the examples

```
cd build/src
./rtty_rx_from_file 8000 < ../../examples/rtty.raw
./rtty_rx_from_file 8000 < ../../examples/rtty.wav
```

To decode a RTTY RAW file in .wav format, you can use 'sox' to convert it first, as follows:

```
sox -t raw -r 8000 -b 16 -e signed -c 1 raw.raw raw.wav
```

To decode the inverse, WAV to RAW:

```
sox rtty.wav --bits 16 rtty.raw
```

# Decoder specs
This decoder expects the two tones to be centered on 2000Hz.
It only decodes 50 baud +/- 225 Hz with stop bit length of 1.5 and no parity bits.

So far tested, with the German DDK agency, decoding works with inverse = True.

More information on stations available and frequencies: https://www.dxinfocentre.com/ratt.htm
Recommended use of kiwiSDR to get signal online.

## Credits

- Dave Freese, W1HKJ for creating fldigi
- Rik van Riel, AB1KW for the RTTY decoder
- Franco Venturi for creating the base of this code for Navtex from which it is based.


## Copyright

(C) 2023 Roger Cano - Licensed under the GNU GPL V3 (see [LICENSE](LICENSE))
