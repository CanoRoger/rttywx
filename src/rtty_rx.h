/*
 * Coyright 2023 Roger Cano
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Code adapted to be used as a library, specially later GNURadio usage.
 * Reduced to receive German RTTY marine forecasts (50 bauds, +/-225Hz
 * shift, 1.5 stop bits, no parity bits) at sea.

// ----------------------------------------------------------------------------
// rtty.h  --  RTTY modem
//
// Copyright (C) 2012
//		Dave Freese, W1HKJ
//		Stefan Fendt, DL1SMF
//
// This file is part of fldigi.
//
// This code bears some resemblance to code contained in gmfsk from which
// it originated.  Much has been changed, but credit should still be
// given to Tomi Manninen (oh2bns@sral.fi), who so graciously distributed
// his gmfsk modem under the GPL.
//
// Fldigi is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Fldigi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fldigi.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------------
*/

#ifndef _RTTY_H
#define _RTTY_H

#include <iostream>

#include "complex.h"
#include "fftfilt.h"
#include <cstdio>
#include <string>
#include <vector>

#define	RTTY_SampleRate	8000
//#define RTTY_SampleRate 11025
//#define RTTY_SampleRate 12000

#define MAXPIPE			1024
#define MAXBITS			(2 * RTTY_SampleRate / 23 + 1)

#define	LETTERS	0x100
#define	FIGURES	0x200c

class Cmovavg {
    //=====================================================================
    // Moving average filter
    //=====================================================================

#define MAXMOVAVG 2048
private:
    double	*in;
    double	out;
    int		len, pint;
    bool	empty;
public:
    Cmovavg(int filtlen = 64);
    ~Cmovavg();
    double run(double a);
    void setLength(int filtlen);
    void reset();
    double value() { return out / (len > 0 ? len : 1); }
};

class fftfilt;

class rtty_rx {
public:
    rtty_rx(int sample_rate, bool only_sitor_b, bool reverse,
            FILE * rawfile=stdout, FILE * messagesfile=nullptr,
            FILE * logfile=stderr);

    int rx_process(const short *buf, int len);

    enum RTTY_RX_STATE {
        RTTY_RX_STATE_IDLE = 0,
        RTTY_RX_STATE_START,
        RTTY_RX_STATE_DATA,
        RTTY_RX_STATE_PARITY,
        RTTY_RX_STATE_STOP,
        RTTY_RX_STATE_STOP2
    };
    enum RTTY_PARITY {
        RTTY_PARITY_NONE = 0,
        RTTY_PARITY_EVEN,
        RTTY_PARITY_ODD,
        RTTY_PARITY_ZERO,
        RTTY_PARITY_ONE
    };

private:
    int m_sample_rate;
    bool m_only_sitor_b;
    bool m_reverse;
    FILE * m_rawfile;
    FILE * m_messagesfile;

    double m_center_frequency_f;

    double shift;
    int symbollen;
    int nbits;
    int stoplen;

    double		rtty_squelch;
    double		rtty_shift;
    double		rtty_BW;
    double		rtty_baud;
    int 		rtty_bits;
    RTTY_PARITY	rtty_parity;
    int			rtty_stop;

    double		mark_noise;
    double		space_noise;
    Cmovavg		*bits;
    bool		nubit;
    bool		bit;

    bool		bit_buf[MAXBITS];

    double mark_phase;
    double space_phase;
    fftfilt *mark_filt;
    fftfilt *space_filt;
    int filter_length;

    double *pipe;
    double *dsppipe;
    int pipeptr;

    cmplx mark_history[MAXPIPE];
    cmplx space_history[MAXPIPE];

    RTTY_RX_STATE rxstate;

    int counter;
    int bitcntr;
    int rxdata;

    double mark_mag;
    double space_mag;
    double mark_env;
    double space_env;
    double	noise_floor;

    unsigned char lastchar;

    int rxmode;
    int shift_state;
    bool rx_lowercase;

    inline cmplx mixer(double &phase, double f, cmplx in);
    unsigned char Bit_reverse(unsigned char in, int n);
    int decode_char();
    bool rx(bool bit);
    char baudot_dec(unsigned char data);
    bool is_mark_space(int &);
    bool is_mark();
    void put_rx_char(int c);
};

#endif
