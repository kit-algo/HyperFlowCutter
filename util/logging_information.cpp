#include "../definitions.h"

namespace LoggingInformation {
	static int g_output_detail = 0;
	void set_output_detail(int __output_detail) { g_output_detail = __output_detail; }
	int get_output_detail() { return g_output_detail; }
}
