/*
    This file is part of spixel.

    Spixel is free software : you can redistribute it and / or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar. If not, see <http://www.gnu.org/licenses/>.
*/


// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#ifdef WIN32
#include "targetver.h"
#endif

#include <stdio.h>
#ifdef WIN32
#include <tchar.h>
#endif
#include <ctime>
#include <iostream>


// TODO: reference additional headers your program requires here
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/photo/photo.hpp>

//[ For finding Memory leaks
// #define _CRTDBG_MAP_ALLOC
// #include <stdlib.h>
// #include <crtdbg.h>
//]
