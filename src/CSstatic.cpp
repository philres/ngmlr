/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Contact: philipp.rescheneder@univie.ac.at
 */

#include "CS.h"


uint CS::prefixBasecount = 13;
uint CS::prefixBits = prefixBasecount * 2;
ulong CS::prefixMask = ((ulong) 1 << prefixBits) - 1;

inline int min(int a, int b) {
	return (a < b) ? a : b;
}

// A->0 C->1 T->2 G->3
inline char encode(char c) {
	return (c >> 1) & 3;
}


// Iteriert ueber jeden Praefix in sequence und fuehrt fuer diesen die Funktion func aus
void CS::PrefixIteration(char const * sequence, uloc length, PrefixIterationFn func, ulong mutateFrom, ulong mutateTo, void* data, uint prefixskip, uloc offset) {
	if (length < prefixBasecount)
		return;

	if (*sequence == 'N') {
		uint n_skip = 1;
		while (*(sequence + n_skip) == 'N')
			++n_skip;

		sequence += n_skip;

		if (n_skip >= (length - prefixBasecount))
			return;
		length -= n_skip;
		offset += n_skip;
	}

	ulong prefix = 0;
	for (uloc i = 0; i < prefixBasecount - 1; ++i) {
		char c = *(sequence + i);
		if (c == 'N') {
			PrefixIteration(sequence + i + 1, length - i - 1, func, mutateFrom, mutateTo, data, prefixskip, offset + i + 1);
			return;
		}

		prefix = prefix << 2;
		char cx = encode(c);
		prefix |= cx;
	}

	uint skipcount = prefixskip;
	for (uloc i = prefixBasecount - 1; i < length; ++i) {
		char c = *(sequence + i);
		if (c == 'N') {
			PrefixIteration(sequence + i + 1, length - i - 1, func, mutateFrom, mutateTo, data, prefixskip, offset + i + 1);
			return;
		}

		prefix = prefix << 2;
		char cx = encode(*(sequence + i));
		prefix |= cx;
		prefix &= prefixMask;

		if (skipcount == prefixskip) {
			func(prefix, offset + i + 1 - prefixBasecount, mutateFrom, mutateTo, data);
			skipcount = 0;
		} else {
			++skipcount;
		}
	}
}
