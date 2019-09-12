#ifndef __ROBOIOSTREAM_HPP__
#define __ROBOIOSTREAM_HPP__

#include <iostream>

class RoboIOStream : public iostream
{
    RoboIOStream(string filename, ios_base::openmode operation, FileTypeFormat type);
};

#endif // __ROBOIOSTREAM_HPP__
