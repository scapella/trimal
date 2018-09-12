//CMake generated code
#include "ReadWriteMS/clustal_state.h"
#include "ReadWriteMS/fasta_m10_state.h"
#include "ReadWriteMS/fasta_state.h"
#include "ReadWriteMS/htmlreport_state.h"
#include "ReadWriteMS/mega_interleaved_state.h"
#include "ReadWriteMS/mega_sequential_state.h"
#include "ReadWriteMS/nexus_m10_state.h"
#include "ReadWriteMS/nexus_state.h"
#include "ReadWriteMS/phylip32_m10_state.h"
#include "ReadWriteMS/phylip32_state.h"
#include "ReadWriteMS/phylip40_m10_state.h"
#include "ReadWriteMS/phylip40_state.h"
#include "ReadWriteMS/phylip_paml_m10_state.h"
#include "ReadWriteMS/phylip_paml_state.h"
#include "ReadWriteMS/pir_state.h"

#include "ReadWriteMS/ReadWriteMachineState.h"

ReadWriteMS::ReadWriteMS()
{
	addState(new clustal_state(this));
	addState(new fasta_m10_state(this));
	addState(new fasta_state(this));
	addState(new htmlreport_state(this));
	addState(new mega_interleaved_state(this));
	addState(new mega_sequential_state(this));
	addState(new nexus_m10_state(this));
	addState(new nexus_state(this));
	addState(new phylip32_m10_state(this));
	addState(new phylip32_state(this));
	addState(new phylip40_m10_state(this));
	addState(new phylip40_state(this));
	addState(new phylip_paml_m10_state(this));
	addState(new phylip_paml_state(this));
	addState(new pir_state(this));
};