/* 046267 Computer Architecture - Spring 21 - HW #3               */
/* Implementation (skeleton)  for the dataflow statistics calculator */

#include <iostream>
#include <vector>
#include <string>
#include "dflow_calc.h"

using namespace std;

#define NUMBER_OF_REGISTERS (32)
#define MAX_OPS (32)

const int DST = (0);
const int SRC1 = (1);
const int SRC2 = (2);
const int SUCCESS = (0);
const int FAILURE = (-1);

int ENTRY_ORDINAL_NUMBER = -1;
int EXIT_ORDINAL_NUMBER = -2;

struct Node {
    int instruction_ordinal_number = ENTRY_ORDINAL_NUMBER;
    unsigned int opcode = -1;
    int instruction_latency = 0;
    int dst_index = -1;
    unsigned int src1_index = ENTRY_ORDINAL_NUMBER;
    unsigned int src2_index = ENTRY_ORDINAL_NUMBER;
    int depth = 0;

    void _add_information(int ordinal_number, const InstInfo instruction_info, const unsigned int opsLatency[]) {
        instruction_ordinal_number = ordinal_number;
        opcode = instruction_info.opcode;
        instruction_latency = (opcode != -1) ? (int)opsLatency[opcode] : 0;
        dst_index = instruction_info.dstIdx;
        src1_index = instruction_info.src1Idx;
        src2_index = instruction_info.src2Idx;
    }
};

typedef vector<vector<Node>> InstructionMatrix;

struct Register {
    bool is_used = false;
    int master_inst_ordinal_number = ENTRY_ORDINAL_NUMBER;

    void _mark_used(int inst_ordinal_num) {
        is_used = true;
        master_inst_ordinal_number = inst_ordinal_num;
    }
};

struct DataFlowGraph {

    InstructionMatrix instruction_vertices;

    // this array will serve to identify the instructions upon which Exit nodes depends. does not include Entry or Exit
    vector<int> pointed_instructions;
    Node *exit_dependencies_array;

    // NOTE: the index in this array represents the OPCODE and not the ORDINAL NUMBER OF THE INSTRUCTION in the trace.
    unsigned int *opcodes_latencies;
    vector<Register> used_registers;
    unsigned int eff_number_of_instructions;

    void _initialize_pointed_instr() {
        for (int i = 0; i < eff_number_of_instructions - 2; ++i) {
            pointed_instructions[i] = 0;
        }
    }

    void _update_latencies(const unsigned int opsLatency[]) {
        opcodes_latencies = new unsigned int[MAX_OPS];

        for (int i = 0; i < MAX_OPS; ++i) {
            opcodes_latencies[i] = opsLatency[i];
        }
    }


    void _init_data_flow_graph(const unsigned int opsLatency[]) {
        // allocate nodes inside and initialise them

        _initialize_pointed_instr();

        _update_latencies(opsLatency);
    }

    DataFlowGraph(const unsigned int opsLatency[], const InstInfo progTrace[], unsigned int numOfInsts) : eff_number_of_instructions(numOfInsts + 2){

        instruction_vertices.resize(eff_number_of_instructions, vector<Node>(3));
        pointed_instructions.resize(numOfInsts);
        used_registers.resize(NUMBER_OF_REGISTERS);
        exit_dependencies_array = new Node[numOfInsts];

        ENTRY_ORDINAL_NUMBER = numOfInsts;
        EXIT_ORDINAL_NUMBER = numOfInsts + 1;

        _init_data_flow_graph(opsLatency);

        _create_dependencies_graph(opsLatency, progTrace);
    }

    ~DataFlowGraph() {
        delete[] exit_dependencies_array;
        delete[] opcodes_latencies;
    }

    void _reset_fields(InstInfo *inst_info) {
        /*** Helper Function to reset all instruction info struct fields. used in _add_dependencies
         *   in case the instruction depends only on Entry node. ***/
        inst_info->opcode = -1;
        inst_info->dstIdx = -1;
        inst_info->src1Idx = ENTRY_ORDINAL_NUMBER;
        inst_info->src2Idx = ENTRY_ORDINAL_NUMBER;
    }

    void _get_instruction_info_by_id(int id, InstInfo *info) {
        info->opcode = instruction_vertices[id][DST].opcode;
        info->dstIdx = instruction_vertices[id][DST].dst_index;
        info->src1Idx = instruction_vertices[id][DST].src1_index;
        info->src2Idx = instruction_vertices[id][DST].src2_index;
    }

    void _add_exit_node_dependencies() {
        /*** Function to handle the dependencies of an exit node. this case is special since we need to find all
         *   the nodes which are not pointed by anyone, i.e., they are the sources of the graph. we need to take
         *   into account the fact the the Entry node is itself a source of the graph, but we can identify it by
         *   it's ordinal number and ignore it. use the pointed instructions array to find which instr. Exit should
         *   point at. ***/

        int counter = 0;
        for (int id = 0; id < eff_number_of_instructions - 2; ++id) {
            if(pointed_instructions[id] == 0) {
                // no one points to this instruction
                InstInfo info;
                _get_instruction_info_by_id(id, &info);
                exit_dependencies_array[counter]._add_information(id,info, opcodes_latencies);
                counter++;
            }
        }
    }

    void _add_dependencies(int inst_ordinal_number, const InstInfo progTrace[]) {

        /*** update sources dependency. if srcX was not written by previous instruction (i.e., no RAW)
             need to point to Entry node. if we're an exit node we need to point to all un-pointed instructions ***/

        InstInfo src1_inst_info, src2_inst_info;
        int src1_ordinal_number = ENTRY_ORDINAL_NUMBER, src2_ordinal_number = ENTRY_ORDINAL_NUMBER;
        // in case there's no RAW, and srcX wasn't written to by previous instruction, point to Entry node.
        _reset_fields(&src1_inst_info);
        _reset_fields(&src2_inst_info);

        // check if there's a dependency and upon which inst. for each one of the srcs
        if(used_registers[progTrace[inst_ordinal_number].src1Idx].is_used) {
            src1_ordinal_number = used_registers[progTrace[inst_ordinal_number].src1Idx].master_inst_ordinal_number;

            src1_inst_info = progTrace[src1_ordinal_number];
            pointed_instructions[src1_ordinal_number] = 1;

        }
        if(used_registers[progTrace[inst_ordinal_number].src2Idx].is_used) {
            src2_ordinal_number = used_registers[progTrace[inst_ordinal_number].src2Idx].master_inst_ordinal_number;

            src2_inst_info = progTrace[src2_ordinal_number];
            pointed_instructions[src2_ordinal_number] = 1;
        }

        // update node's dependencies.
        instruction_vertices[inst_ordinal_number][SRC1]._add_information(src1_ordinal_number, src1_inst_info, opcodes_latencies);
        instruction_vertices[inst_ordinal_number][SRC2]._add_information(src2_ordinal_number, src2_inst_info, opcodes_latencies);
    }

    int _get_latency(int i) {
        return instruction_vertices[i][DST].depth + instruction_vertices[i][DST].instruction_latency;
    }

    int _calc_depth(int instruction_ordinal_number) {
        if (instruction_ordinal_number == ENTRY_ORDINAL_NUMBER)
            return 0;

        int dep_a_depth = _get_latency(instruction_vertices[instruction_ordinal_number][SRC1].instruction_ordinal_number);
        int dep_b_depth = _get_latency(instruction_vertices[instruction_ordinal_number][SRC2].instruction_ordinal_number);
        instruction_vertices[instruction_ordinal_number][DST].depth = max(dep_a_depth, dep_b_depth);
        return instruction_vertices[instruction_ordinal_number][DST].depth;
    }

    int _update_intrcutions_depths() {
        for (int i = 0; i < eff_number_of_instructions - 2; ++i) {
            _calc_depth(i);
        }
    }

    void _add_vertex(int instruction_ordinal_number, const InstInfo progTrace[]) {
        /*** Function to update specific instruction by it's ordinal number.
         *   for the i'th instruction we will update the adequate node in the instruction_vertices vector
         *   and update the adjacency list holding it's dependencies. ***/

        instruction_vertices[instruction_ordinal_number][DST]._add_information(instruction_ordinal_number, progTrace[instruction_ordinal_number], opcodes_latencies);

        _add_dependencies(instruction_ordinal_number, progTrace);
    }

    void _create_dependencies_graph(const unsigned int opsLatency[], const InstInfo progTrace[]) {
        /*** Function to create the dependencies graph. loop over all instuctions and in ordinal order do:
         *      1. add vertex and dependencies
         *      2. mark dst register as used ***/

        int used_register_number;

        for (int inst_ordinal_num = 0; inst_ordinal_num < eff_number_of_instructions - 2; ++inst_ordinal_num) {

            // assume this function will also add the dependencies
            _add_vertex(inst_ordinal_num, progTrace);

            used_register_number = progTrace[inst_ordinal_num].dstIdx;
            used_registers[used_register_number]._mark_used(inst_ordinal_num);
        }

        _update_intrcutions_depths();

        _add_exit_node_dependencies();

        instruction_vertices[ENTRY_ORDINAL_NUMBER][DST].instruction_ordinal_number = ENTRY_ORDINAL_NUMBER;
        instruction_vertices[EXIT_ORDINAL_NUMBER][DST].instruction_ordinal_number = EXIT_ORDINAL_NUMBER;
    }


    int getProgDepthRec(){
        int idx = 0, max_depth = 0, curr_depth = 0, latency = 0;
        while (exit_dependencies_array[idx].instruction_ordinal_number != -1) {
            curr_depth = instruction_vertices[exit_dependencies_array[idx].instruction_ordinal_number][DST].depth +
                         instruction_vertices[exit_dependencies_array[idx].instruction_ordinal_number][DST].instruction_latency;
            if(curr_depth > max_depth) {
                max_depth = curr_depth;
            }
            idx++;
        }
        return max_depth;
    }


    int getInstDeps(unsigned int theInst, int *src1DepInst, int *src2DepInst) {
        *src1DepInst = instruction_vertices[theInst][SRC1].instruction_ordinal_number;
        if(*src1DepInst == ENTRY_ORDINAL_NUMBER) *src1DepInst = -1;
        *src2DepInst = instruction_vertices[theInst][SRC2].instruction_ordinal_number;
        if(*src2DepInst == ENTRY_ORDINAL_NUMBER) *src2DepInst = -1;
        return SUCCESS;
    }
};


ProgCtx analyzeProg(const unsigned int opsLatency[], const InstInfo progTrace[], unsigned int numOfInsts) {
    try {
        DataFlowGraph *ctx = new DataFlowGraph(opsLatency, progTrace, numOfInsts);
        return ctx;
    }

    catch(std::bad_alloc& e){
        return PROG_CTX_NULL;
    }
}

void freeProgCtx(ProgCtx ctx) {
    if(!ctx) return;
    delete (DataFlowGraph*)ctx;
}

int getInstDepth(ProgCtx ctx, unsigned int theInst) {
    if(ctx == nullptr) return -1;

    int latency = 0;
    DataFlowGraph *dflow = (DataFlowGraph*)ctx;

    if (theInst < 0 || theInst >= (dflow->eff_number_of_instructions-2)) {
        return FAILURE;
    }
    return dflow->instruction_vertices[theInst][DST].depth;
}

int getInstDeps(ProgCtx ctx, unsigned int theInst, int *src1DepInst, int *src2DepInst) {
    if (!ctx || !src1DepInst || !src2DepInst) return FAILURE;
    DataFlowGraph *dflow = (DataFlowGraph*)ctx;

    if (theInst < 0 || theInst >= (dflow->eff_number_of_instructions-2)) {
        return FAILURE;
    }

    return dflow->getInstDeps(theInst, src1DepInst, src2DepInst);
}

int getProgDepth(ProgCtx ctx) {
    if(ctx == nullptr) return -1;
    DataFlowGraph *dflow = (DataFlowGraph*)ctx;

    // pass the Exit node ordinal number as the source ordinal number from which to start calculating the path.
    return dflow->getProgDepthRec();
}


