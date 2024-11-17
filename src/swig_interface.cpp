#include "swig_interface.h"
#include "contextsv.h"

/// @cond
#include <iostream>
/// @endcond


// Run the CLI with the given parameters
int run(InputData input_data)
{

	// Run ContextSV
	ContextSV contextsv(input_data);
	try
	{	
		contextsv.run();
	}

	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return -1;
	}

	return 0;
}
