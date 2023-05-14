/*
 * Copyright 2023 Roger Cano
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Code adapted to be used as a library, specially later GNURadio usage.
 * Reduced to receive German RTTY marine forecasts (50 bauds, +/-225Hz
 * shift, 1.5 stop bits, no parity bits) at sea.
 */

/*
// ----------------------------------------------------------------------------
// rtty.cxx  --  RTTY modem
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

#include <iostream>
#include <fstream>
//#include <stdio.h>
#include <climits>
#include <cstring>

#include "misc.h"
#include "fftfilt.h"
#include "threads.h"
#include "rtty_rx.h"


#define FILTER_DEBUG 0
#define SHAPER_BAUD 150

#define TWOPI 3.14169

static const int deviation_f = 225; // half frequency shift of the rtty
static const double dflt_center_freq = 2000.0 ; // central frequency of the rtty

// Minimum length of logged messages
static const size_t min_siz_logged_msg = 0;

static char letters[32] = {
    '\0',	'E',	'\n',	'A',	' ',	'S',	'I',	'U',
    '\r',	'D',	'R',	'J',	'N',	'F',	'C',	'K',
    'T',	'Z',	'L',	'W',	'H',	'Y',	'P',	'Q',
    'O',	'B',	'G',	' ',	'M',	'X',	'V',	' '
};
static char figures[32] = {
    '\0',	'3',	'\n',	'-',	' ',	'\a',	'8',	'7',
    '\r',	'$',	'4',	'\'',	',',	'!',	':',	'(',
    '5',	'"',	')',	'2',	'#',	'6',	'0',	'1',
    '9',	'?',	'&',	' ',	'.',	'/',	';',	' '
};

int dspcnt = 0;

static char msg1[20];

// printing the result
void rtty_rx::put_rx_char(int c) {
    if (m_rawfile != nullptr)
    putc(c, m_rawfile);
}

rtty_rx::rtty_rx(int sample_rate, bool only_sitor_b, bool reverse,
                 FILE * rawfile, FILE * messagesfile, FILE * logfile) {
    //printf("rtty_rx-start 2");
    m_sample_rate = sample_rate;
    m_only_sitor_b = only_sitor_b;
    m_reverse = reverse;
    m_rawfile = rawfile;
    m_messagesfile = messagesfile;

    m_center_frequency_f = dflt_center_freq;
    rx_lowercase = 0;
    rtty_shift = 450;
    shift = 450;
    rtty_baud = 50;
    filter_length = 512;

    // old rx_init():

    rxstate = RTTY_RX_STATE_IDLE;

//    for (int i = 0; i < MAXBITS; i++ ) bit_buf[i] = 0.0;

    mark_phase = 0;
    space_phase = 0;

    mark_mag = 0;
    space_mag = 0;
    mark_env = 0;
    space_env = 0;

//    inp_ptr = 0;
    lastchar = 0;

    // if (mark_filt) delete mark_filt; // these caused segmentation fault, exit -11 in gnuradio
    mark_filt = new fftfilt(rtty_baud/m_sample_rate, filter_length);
    mark_filt->rtty_filter(rtty_baud/m_sample_rate);

    //if (space_filt) delete space_filt;
    space_filt = new fftfilt(rtty_baud/m_sample_rate, filter_length);
    space_filt->rtty_filter(rtty_baud/m_sample_rate);

    for (int i = 0; i < MAXPIPE; i++) mark_history[i] = space_history[i] = cmplx(0,0);

    lastchar = 0;
    // Copied from old rtty_rx()    printf("halfway of rx_start");


    bits = (Cmovavg *)0;

    pipe = new double[MAXPIPE];
    dsppipe = new double [MAXPIPE];

    // Copied from old reset():

    double stl;

    rtty_bits = 5;
    nbits = rtty_bits;
    rtty_parity = RTTY_PARITY_NONE;

    shift_state = LETTERS;
    rxmode = LETTERS;

    symbollen = (int) (m_sample_rate / rtty_baud + 0.5);

    rtty_BW = rtty_baud * 2;

//  Copied from reset_filters();

    if (bits)
        bits->setLength(symbollen / 8);//2);
    else
        bits = new Cmovavg(symbollen / 8);//2);

    mark_noise = space_noise = 0;
    bit = nubit = true;

// Stop symbol length 1.5 for marine RTTY weather forecast
    stl = 1.5;

    stoplen = (int) (stl * sample_rate / rtty_baud + 0.5);
    pipeptr = 0;

//    for (int i = 0; i < MAXBITS; i++ ) bit_buf[i] = 0.0;

    dspcnt = 2*(nbits + 2);

//    clear_zdata = true;
    for (int i = 0; i < MAXPIPE; i++) mark_history[i] = space_history[i] = cmplx(0,0);
    //printf("end of rx_start");
}

cmplx rtty_rx::mixer(double &phase, double f, cmplx in)
{
    cmplx z = cmplx( cos(phase), sin(phase)) * in;

    phase -= TWOPI * f / m_sample_rate; //sample_rate;
    if (phase < -TWOPI) phase += TWOPI;

    return z;
}

Cmovavg::Cmovavg (int filtlen)
{
    //=====================================================================
    // Moving average filter
    //
    // Simple in concept, sublime in implementation ... the fastest filter
    // in the west.  Also optimal for the processing of time domain signals
    // characterized by a transition edge.  The is the perfect signal filter
    // for CW, RTTY and other signals of that type.  For a given filter size
    // it provides the greatest s/n improvement while retaining the sharpest
    // leading edge on the filtered signal.
    //=====================================================================
    len = filtlen;
    in = new double[len];
    empty = true;
}

Cmovavg::~Cmovavg()
{
    if (in) delete [] in;
}

double Cmovavg::run(double a)
{
    if (!in) {
        return a;
    }
    if (empty) {
        empty = false;
        out = 0;
        for (int i = 0; i < len; i++) {
            in[i] = a;
            out += a;
        }
        pint = 0;
        return a;
    }
    out = out - in[pint] + a;
    in[pint] = a;
    if (++pint >= len) pint = 0;
    return out / len;
}

void Cmovavg::setLength(int filtlen)
{
    if (filtlen > len) {
        if (in) delete [] in;
        in = new double[filtlen];
    }
    len = filtlen;
    empty = true;
}

void Cmovavg::reset()
{
    empty = true;
}

unsigned char rtty_rx::Bit_reverse(unsigned char in, int n)
{
    unsigned char out = 0;

    for (int i = 0; i < n; i++)
        out = (out << 1) | ((in >> i) & 1);

    return out;
}

static int rparity(int c)
{
    int w = c;
    int p = 0;
    while (w) {
        p += (w & 1);
        w >>= 1;
    }
    return p & 1;
}

int rttyparity(unsigned int c, int nbits)
{
    c &= (1 << nbits) - 1;


    // No parity for weather RTTY

//    switch (progdefaults.rtty_parity) {
//    default:
//    case rtty_rx::RTTY_PARITY_NONE:
        return 0;

//    case rtty_rx::RTTY_PARITY_ODD:
//        return rparity(c);

//    case rtty_rx::RTTY_PARITY_EVEN:
//        return !rparity(c);

//    case rtty_rx::RTTY_PARITY_ZERO:
//        return 0;

//    case rtty_rx::RTTY_PARITY_ONE:
//        return 1;
//    }
}

int rtty_rx::decode_char()
{
    unsigned int parbit, par, data;

    parbit = (rxdata >> nbits) & 1;
    par = rttyparity(rxdata, nbits);

    if (rtty_parity != RTTY_PARITY_NONE && parbit != par)
        return 0;

    data = rxdata & ((1 << nbits) - 1);

    if (nbits == 5)
        return baudot_dec(data);

    return data;
}

bool rtty_rx::is_mark_space( int &correction)
{
    correction = 0;

// test for rough bit position
    if (bit_buf[0] && !bit_buf[symbollen-1]) {

// test for mark/space straddle point
        for (int i = 0; i < symbollen; i++)
            correction += bit_buf[i];

        if (abs(symbollen/2 - correction) < 6) // too small & bad signals are not decoded
            return true;
    }
    return false;
}

bool rtty_rx::is_mark()
{
    return bit_buf[symbollen / 2];
}

bool rtty_rx::rx(bool bit) // original modified for probability test
{
    bool flag = false;
    unsigned char c = 0;
    int correction;

    for (int i = 1; i < symbollen; i++) bit_buf[i-1] = bit_buf[i];
    bit_buf[symbollen - 1] = bit;

    switch (rxstate) {
    case RTTY_RX_STATE_IDLE:
        if ( is_mark_space(correction)) {
            rxstate = RTTY_RX_STATE_START;
            counter = correction;
        }
        break;
    case RTTY_RX_STATE_START:
        if (--counter == 0) {
            if (!is_mark()) {
                rxstate = RTTY_RX_STATE_DATA;
                counter = symbollen;
                bitcntr = 0;
                rxdata = 0;
            } else {
                rxstate = RTTY_RX_STATE_IDLE;
            }
        }
        break;
    case RTTY_RX_STATE_DATA:
        if (--counter == 0) {
            rxdata |= is_mark() << bitcntr++;
            counter = symbollen;
        }
        if (bitcntr == nbits + (rtty_parity != RTTY_PARITY_NONE ? 1 : 0))
            rxstate = RTTY_RX_STATE_STOP;
        break;
    case RTTY_RX_STATE_STOP:
        if (--counter == 0) {
            if (is_mark()) {
                if (true) {
                    c = decode_char();
                    if ( c != 0 ) {
                        // supress <CR><CR> and <LF><LF> sequences
                        // these were observed during the RTTY contest 2/9/2013
                        if (c == '\r' && lastchar == '\r');
                        else if (c == '\n' && lastchar == '\n');
                        else
                            // put_rx_char gives the final data - end of processing
                            put_rx_char(rx_lowercase ? tolower(c) : c);

                        lastchar = c;
                    }
                    flag = true;
                }
            }
            rxstate = RTTY_RX_STATE_IDLE;
        }
        break;
    default : break;
    }
    return flag;
}

int rtty_rx::rx_process(const short *buf, int len)
{
    //printf("in rx_process");
    const short *buffer = buf;
    int length = len;
//    static int showxy = symbollen;

    cmplx z, zmark, zspace, *zp_mark, *zp_space;

    int n_out = 0;
    static int bitcount = 5 * nbits * symbollen;

    while (length-- > 0) {

// Create analytic signal from sound card input samples

    z = cmplx(*buffer, *buffer);
    buffer++;

// Mix it with the audio carrier frequency to create two baseband signals
// mark and space are separated and processed independently
// lowpass Windowed Sinc - Overlap-Add convolution filters.
// The two fftfilt's are the same size and processed in sync
// therefore the mark and space filters will concurrently have the
// same
    //#endifsize outputs available for further processing

        zmark = mixer(mark_phase, m_center_frequency_f + deviation_f, z);
        mark_filt->run(zmark, &zp_mark);

        zspace = mixer(space_phase, m_center_frequency_f - deviation_f, z);
        n_out = space_filt->run(zspace, &zp_space);

// envelope & noise levels for mark & space
        for (int i = 0; i < n_out; i++) {

// determine noise floor & envelope for mark & space
            mark_mag = abs(zp_mark[i]);
            mark_env = decayavg (mark_env, mark_mag,
                        (mark_mag > mark_env) ? symbollen / 4 : symbollen * 16);
            mark_noise = decayavg (mark_noise, mark_mag,
                        (mark_mag < mark_noise) ? symbollen / 4 : symbollen * 48);
            space_mag = abs(zp_space[i]);
            space_env = decayavg (space_env, space_mag,
                        (space_mag > space_env) ? symbollen / 4 : symbollen * 16);
            space_noise = decayavg (space_noise, space_mag,
                        (space_mag < space_noise) ? symbollen / 4 : symbollen * 48);


            noise_floor = std::min(space_noise, mark_noise);

// clipped if clipped decoder selected
// mclipped aka. mark_abs // sclipped aka. space_abs

            double mclipped = 0, sclipped = 0;
            mclipped = mark_mag > mark_env ? mark_env : mark_mag;
            sclipped = space_mag > space_env ? space_env : space_mag;
            if (mclipped < noise_floor) mclipped = noise_floor;
            if (sclipped < noise_floor) sclipped = noise_floor;


// mark-space discriminator with automatic threshold
// correction, see:
// http://www.w7ay.net/site/Technical/ATC/

//			double v0, v1, v2, v3, v4, v5;
            double v3;

// no ATC
//            v0 = mark_mag - space_mag;
// Linear ATC
//			v1 = mark_mag - space_mag - 0.5 * (mark_env - space_env);
// Clipped ATC
//			v2  = (mclipped - noise_floor) - (sclipped - noise_floor) - 0.5 * (
//					(mark_env - noise_floor) - (space_env - noise_floor));
// Optimal ATC
            v3  = (mclipped - noise_floor) * (mark_env - noise_floor) -
                    (sclipped - noise_floor) * (space_env - noise_floor) - 0.25 * (
                    (mark_env - noise_floor) * (mark_env - noise_floor) -
                    (space_env - noise_floor) * (space_env - noise_floor));
// Kahn Squarer with Linear ATC
//			v4 =  (mark_mag - noise_floor) * (mark_mag - noise_floor) -
//					(space_mag - noise_floor) * (space_mag - noise_floor) - 0.25 * (
//					(mark_env - noise_floor) * (mark_env - noise_floor) -
//					(space_env - noise_floor) * (space_env - noise_floor));
// Kahn Squarer with Clipped ATC
//			v5 =  (mclipped - noise_floor) * (mclipped - noise_floor) -
//					(sclipped - noise_floor) * (sclipped - noise_floor) - 0.25 * (
//					(mark_env - noise_floor) * (mark_env - noise_floor) -
//					(space_env - noise_floor) * (space_env - noise_floor));
//				switch (progdefaults.rtty_demodulator) {
//			switch (2) { // Optimal ATC
//			case 0: // linear ATC
//				bit = v1 > 0;
//				break;
//			case 1: // clipped ATC
//				bit = v2 > 0;
//				break;
//			case 2: // optimal ATC
                bit = v3 > 0;
//				break;
//			case 3: // Kahn linear ATC
//				bit = v4 > 0;
//				break;
//			case 4: // Kahn clipped
//				bit = v5 > 0;
//				break;
//			case 5: // No ATC
//			default :
//                bit = v0 > 0;
//			}

// detect TTY signal transitions
// rx(...) returns true if valid TTY bit stream detected
// either character or idle signal

                if ( rx( m_reverse ? !bit : bit ) ) {        //rx( reverse ? !bit : bit ) ) {
                    dspcnt = symbollen * (nbits + 2);

//                    clear_zdata = true;
                    bitcount = 5 * nbits * symbollen;

                    if (bitcount) --bitcount;
                }
        }
    }
    return 0;
}

char rtty_rx::baudot_dec(unsigned char data)
{
    int out = 0;

    switch (data) {
    case 0x1F:		/* letters */
        rxmode = LETTERS;
        break;
    case 0x1B:		/* figures */
        rxmode = FIGURES;
        break;
    case 0x04:		/* unshift-on-space */
//        if (progdefaults.UOSrx)
//            rxmode = LETTERS;
        return ' ';
        break;
    default:
        if (rxmode == LETTERS)
            out = letters[data];
        else
            out = figures[data];
        break;
    }

    return out;
}
