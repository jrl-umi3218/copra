/*      File: api.h
*       This file is part of the program copra
*       Program description : This a C++ implementation of a Time Invariant Linear Model Predictive Controller (LMPC) done in C++14 with python bindings
*       Copyright (C) 2017 -  Vincent Samy (LIRMM). All Right reserved.
*
*       This software is free software: you can redistribute it and/or modify
*       it under the terms of the CeCILL-C license as published by
*       the CEA CNRS INRIA, either version 1
*       of the License, or (at your option) any later version.
*       This software is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*       CeCILL-C License for more details.
*
*       You should have received a copy of the CeCILL-C License
*       along with this software. If not, it can be found on the official website
*       of the CeCILL licenses family (http://www.cecill.info/index.en.html).
*/
#pragma once

// Handle portable symbol export.
// Defining manually which symbol should be exported is required
// under Windows whether MinGW or MSVC is used.
//
// The headers then have to be able to work in two different modes:
// - dllexport when one is building the library,
// - dllimport for clients using the library.
//
// On Linux, set the visibility accordingly. If C++ symbol visibility
// is handled by the compiler, see: http://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
// On Microsoft Windows, use dllimport and dllexport to tag symbols.
#define COPRA_DLLIMPORT __declspec(dllimport)
#define COPRA_DLLEXPORT __declspec(dllexport)
#define COPRA_DLLLOCAL
#else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#if __GNUC__ >= 4
#define COPRA_DLLIMPORT __attribute__((visibility("default")))
#define COPRA_DLLEXPORT __attribute__((visibility("default")))
#define COPRA_DLLLOCAL __attribute__((visibility("hidden")))
#else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#define COPRA_DLLIMPORT
#define COPRA_DLLEXPORT
#define COPRA_DLLLOCAL
#endif // __GNUC__ >= 4
#endif // defined _WIN32 || defined __CYGWIN__

#ifdef COPRA_STATIC
// If one is using the library statically, get rid of
// extra information.
#define COPRA_DLLAPI
#define COPRA_LOCAL
#else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#ifdef copra_EXPORTS
#define COPRA_DLLAPI COPRA_DLLEXPORT
#else
#define COPRA_DLLAPI COPRA_DLLIMPORT
#endif // COPRA_EXPORTS
#define COPRA_LOCAL COPRA_DLLLOCAL
#endif // COPRA_STATIC
