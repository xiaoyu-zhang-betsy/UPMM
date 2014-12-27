/*
    This file is part of a demo implementation of an importance sampling technique
    described in the "On-line Learning of Parametric Mixture Models for Light Transport Simulation"
    (SIGGRAPH 2014) paper.
    The implementation is based on Mitsuba, a physically based rendering system.

    Copyright (c) 2014 by Jiri Vorba, Ondrej Karlik, Martin Sik.
    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#if defined( DISCRETE_PHOTON )
/// In the main branch of Mitsuba this is the usual implementation as memory is saved.
#include "discretephoton.cpp"
#else
/// For guiding, however, it is advantageous to incorporate more information for training of the sampling distributions.
#include "nondiscretephoton.cpp"
#endif