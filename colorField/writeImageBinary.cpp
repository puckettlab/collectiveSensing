#include <fstream>
#include <iterator>
#include <algorithm>

const unsigned char * const ones( int size )
{
    unsigned char * data = new unsigned char[size];
    for( int i=0; i<size; ++i )
        data[i] = 255;
    return data;
}


int writeImageBinary()
{
    std::ifstream input( "C:\\Final.gif", std::ios::binary );
    std::ofstream output( "C:\\myfile.gif", std::ios::binary );

    std::copy( 
        std::istreambuf_iterator<char>(input), 
        std::istreambuf_iterator<char>( ),
        std::ostreambuf_iterator<char>(output));
}
