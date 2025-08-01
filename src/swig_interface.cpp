#include "swig_interface.h"
#include "contextsv.h"

/// @cond
#include <iostream>
/// @endcond


// Run the CLI with the given parameters
int run(const InputData& input_data)
{
	// Run ContextSV
	ContextSV contextsv;
	try
	{	
		contextsv.run(input_data);
	}

	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return -1;
	}

	return 0;
}
