//fm_trial.cpp
//Fidducia-Mattheyses method to implement bi-partitioning
//Magicgirl_qian@SUSTech

#include <stdio.h>
#include <stdlib.h>  
#include <stdint.h>
#include <string.h>
#include <string>
#include <list>
#include <map> 
#include <vector> 
#include <iterator>
#include <iostream>
#include <ctime>

using namespace std;

class NET;
class CELL;

//  DECLARE GLOBAL VARIABLES //

multimap <int, CELL*> gain_bucket [2];

vector <CELL*> cells; //pointers to all cells in the design

map <size_t, CELL*> cells_map; //map of cells. Key = cell_id, value = cell pointer
				//useful for accessing cells by name

int cut_size = 0; //total cuts across partitions	
			
size_t partition_area[2] = {0, 0}; //used to track area of partitions

vector <NET*> nets; //pointers to all nets in the design

int min_cut_size;
int base_partition = 0;

bool positive_gains_exhausted = false; //true after the initial gain bucket has been exhausted of positive gains
						//used for updating the best_partition

size_t ratio_factor = 50;
float floor_ratio;


//***************************//
//     FUNCTION HEADERS	     //
//***************************//

void readInput();
void computeGains();
void initialPartition();
void findCutsize();
void saveBestSolution();
void recallBestSolution();
bool FMPartitionPass();
bool performNextMove();
bool cellCanMove(CELL* base_cell);
void printResult();

class NET {
	public:
	size_t net_id;
	size_t num_pins; //total number of pins on this net
    size_t cost;
	bool has_locked[2];//a net is "dead" if no legal cell swaps remain that could change the cutstate
			//this occurs when there is a locked or fixed cell in both partitions on this net
			
	size_t partition_count[2]; //count of how many cells in each partition on this net
	
	vector <CELL*> cell_list; //pointers to all cells on this net			
};

class CELL {
	public:
	size_t cell_id;
	size_t area; 	//physical area of the cell
	bool locked; 	//true for cells locked for the rest of the pass
	int gain; 	//gain is the change in cutset if this cell were to swap partitions
	size_t partition; 	//which partition this cell is in: 0 or 1
	size_t best_partition; //used to remember this cell's partition in the best solution found so far
	multimap <int, CELL*> :: iterator gain_itr; //an iterator that points to this cell's location in the gain buckets		
	
	vector <NET*> net_list; //pointers to all nets this cell is on

	//change the partition of the cell. Update gains for this cell and all neighbors. Update cutsize.
	void changePartition()
	{		
		//remove the cell from the bucket i.e. "lock" the cell until the next pass
		gain_bucket[partition].erase(gain_itr);				
		locked = true;	
		
		//update cut_size. cutsize is reduced by the gain of the moved cell.
		cut_size -= gain;		
			
		int base_partition = partition;
		
		//change partition for this cell
		partition_area[partition] -= area; //decrease area in old partition
		partition = !partition;		
		partition_area[partition] += area; //increase area in new partition
		
		if(gain < 0)
			positive_gains_exhausted = true;
											
		//Update the gain of all neighboring cells		
		//Loop through all nets the cell is on. Update gains for neighboring cells on critical nets
		for(int i = 0; i < net_list.size(); i++)
		{						
			//if the net has a locked (or fixed) cell in both partitions,
			//then the net is "dead" and cannot cause any change in cutstate (prevents wasting time)
			if(net_list[i]->has_locked[0] && net_list[i]->has_locked[1])
			{				
				continue;
			}
							
			//create "From" and "To" partition counts
			size_t from = net_list[i]->partition_count[!partition];   //i.e. F(n), where the cell came from
			size_t to   = net_list[i]->partition_count[partition];    //i.e. T(n), where the cell is going to
			
		//check critical nets before the move
			if(to == 0)
			{
			//increment gain on all free cells on this net
				for(int k = 0; k < net_list[i]->cell_list.size(); k++)
				{
					if((!net_list[i]->cell_list[k]->locked))		
					{
						net_list[i]->cell_list[k]->updateGain(net_list[i]->cost);
					}
				}
			}								
			else if (to == 1)
			{
			//decrement gain of the only "to" cell
				for(int k = 0; k < net_list[i]->cell_list.size(); k++)
				{
					if(!net_list[i]->cell_list[k]->locked && base_partition != net_list[i]->cell_list[k]->partition)
					{						
						net_list[i]->cell_list[k]->updateGain(-1*(int)net_list[i]->cost);
						break;
					}
				}
			}
							
			//update from and to counts
			from--;
			net_list[i]->partition_count[!partition]--;
			to++;
			net_list[i]->partition_count[partition]++;
			
			//check critical nets after the move
			if(from == 0)
			{ 
			//decrement gain on all free cells on this net
				for(int k = 0; k < net_list[i]->cell_list.size(); k++)
				{
					if(!net_list[i]->cell_list[k]->locked)		
					{
						net_list[i]->cell_list[k]->updateGain(-1*(int)net_list[i]->cost);
					}
				}
			}								
			else if (from == 1)
			{
			//increment gain of the only "from" cell
				for(int k = 0; k < net_list[i]->cell_list.size(); k++)
				{
					if(!net_list[i]->cell_list[k]->locked && base_partition == net_list[i]->cell_list[k]->partition)
					{
						net_list[i]->cell_list[k]->updateGain(net_list[i]->cost);
						break;
					}
				}
			}
			
			//update the has_locked member for this net for the "To" partition
			net_list[i]->has_locked[!base_partition] = true;
		}	
	}
	
	
	void updateGain(int gain_change)
	{		
		//if the cell is locked, do nothing. Only update free cells.
		//because locked cells are no longer considered for this pass
		if(locked) return;
		
		//remove the old entry in gain bucket
		gain_bucket[partition].erase(gain_itr);
		
		gain += gain_change;
						
		//add the new entry to gain bucket, and set new gain_itr
		gain_itr = gain_bucket[partition].insert(pair <int, CELL*> (gain, this));		
	}
};


//reads the input of cells and nets
void readInput()
{   
    size_t cell_num;
    cout << "Please enter the number of nodes:";
    cin >> cell_num;

    size_t cell_id, area, max_area = 0, total_area = 0;
    cout << "Please enter each of the " << cell_num << " nodes with its id and the node area:\n";
    for(size_t i = 0; i < cell_num; i++){
        cin >> cell_id >> area;
		if(area > max_area) max_area = area;
		total_area += area;

        CELL* current_cell;
		current_cell = new CELL;

        //vector <NET*> current_net_list; //Is it necessary?
        current_cell->cell_id = cell_id;
        current_cell->area = area;
        current_cell->gain = 0;
        current_cell->locked = false;
        //current_cell->net_list = current_net_list;	
        current_cell->partition = 0;//place all cells in partition 0 initially
        current_cell->best_partition = 0;

        partition_area[0] += current_cell->area;
        cells.push_back(current_cell);
		cells_map[current_cell->cell_id] = current_cell;
    }

    size_t net_num;
    cout << "Please enter the number of edges:";
    cin >> net_num;
    cout << "Please enter each of the "<< net_num <<" edges with the number of connected nodes and their node ids, followed by the edge cost:\n";
    for(size_t i = 0; i < net_num; i++){
        size_t numPins;
        cin >> numPins;
        NET* current_net;
		current_net = new NET;

        current_net->net_id = i;
        current_net->num_pins = numPins;
        for(size_t j = 0; j < numPins; j++){
            cin >> cell_id;
            if(cells_map[cell_id]->net_list.empty() || cells_map[cell_id]->net_list.back() != current_net ){//avoid adding duplicates
                cells_map[cell_id]->net_list.push_back(current_net);
				current_net->cell_list.push_back(cells_map[cell_id]);
            }
        }
        size_t cost;
        cin >> cost;
        current_net->cost = cost;
        nets.push_back(current_net);
    }
    cout << "Please enter the percentage of the ratio factor:";
    cin >> ratio_factor;
	floor_ratio = (float)ratio_factor/100 - (float)max_area/total_area;
}

void initialPartition()
{
	srand(time(0)); //use current time for random seed				
	int index;
	//randomly move cells into other partition until approximately half area in each
	while((float)partition_area[1] / (float)(partition_area[0] + partition_area[1]) < (float)ratio_factor/100 )
	{
		//select a random cell to move to partition 1
		index = rand() % cells.size();
		if(cells[index]->partition == 0)
		{	
			cells[index]->partition = 1;
			cells[index]->best_partition = 1;
			partition_area[0] -= cells[index]->area;
			partition_area[1] += cells[index]->area;		
		}
	}
	cout << "The Initial ";
	findCutsize();
}


//computes the cutsize of the partitions and stores in the global variable "cut_size"
void findCutsize()
{
	cut_size = 0;	
	
	for(size_t i = 0; i < nets.size(); i++)
	{		
		nets[i]->partition_count[0] = 0;
		nets[i]->partition_count[1] = 0;
		
		for(size_t j = 0; j < nets[i]->cell_list.size(); j++)
		{
			nets[i]->partition_count[nets[i]->cell_list[j]->partition]++;			
		}
		
		if(nets[i]->partition_count[0] > 0 && nets[i]->partition_count[1] > 0)
			cut_size += nets[i]->cost;
	}
	cout << "Cutcost:" << cut_size << endl;
}

//runs through all nets in the design. Computes the cell's potential swap gain based on neighbors
//gain is the decrease in cutsize that would be obtained by swapping a cell
void computeGains()
{
	gain_bucket[0].clear();
	gain_bucket[1].clear();
	
	//look through all cells, compute gain for each cell 
	for(size_t i = 0; i < cells.size(); i++)
	{	
		cells[i]->gain = 0;
		
		//for each net the cell is in, check if there is any gain for that net
		for (size_t j = 0; j < cells[i]->net_list.size(); j++)
		{	
			NET* current_net = cells[i]->net_list[j];				
						
			//if there is exactly one cell in the "From" partition, increase gain
			if(current_net->partition_count[cells[i]->partition] == 1)
				cells[i]->gain += current_net->cost;
			
			//if there is exactly zero cells in the "To" partition, decrease gain
			if(current_net->partition_count[!cells[i]->partition] == 0)	
				cells[i]->gain -= current_net->cost;			
		}
		
		//add this cell to a gain bucket, and initialize gain_itr				
		cells[i]->gain_itr = gain_bucket[cells[i]->partition].insert(pair <int, CELL*> (cells[i]->gain,  cells[i]));	
	} 			
}

void saveBestSolution()
{	
	//keep track of the best solution found
	min_cut_size = cut_size;
			
	for(int i = 0; i < cells.size(); i++)
		cells[i]->best_partition = cells[i]->partition;	
}

//run after each pass
//resets the state of the cells to the best solution found so far
//also resets other attributes to set up for the next pass
void recallBestSolution()
{
	cut_size = min_cut_size;
	
	//reset the cells partitions to the best solution
	//unlock all cells
	for(int i = 0; i < cells.size(); i++)
	{
		cells[i]->partition = cells[i]->best_partition;
		cells[i]->locked = false;		
	}				

	//reset the partition counts for the best solution
	for(int i = 0; i < nets.size(); i++)
	{	
		//reset nets	
		nets[i]->partition_count[0] = 0;
		nets[i]->partition_count[1] = 0;
		
		nets[i]->has_locked[0] = 0; 
		nets[i]->has_locked[1] = 0;
		
		for(int j = 0; j < nets[i]->cell_list.size(); j++)
		{
			//count how many cells in each partition on this net
			nets[i]->partition_count[nets[i]->cell_list[j]->partition]++;			
		}		
	}
	
	//recalculate the gains for the best solution
	computeGains(); 
	
	//set the partition areas to the best solution
	partition_area[0] = 0;
	partition_area[1] = 0;
	
	for(int i = 0; i < cells.size(); i++)
	{		
		partition_area[cells[i]->partition] += cells[i]->area;
	}	
}

bool FMPartitionPass()
{
	int starting_cut_size = cut_size;
	positive_gains_exhausted = false; //true after the initial gain bucket has been exhausted of positive gains
					//used to prevent updating the best_partition too often	
	
	//perform moves until no legal moves remain
	while(performNextMove())
	{
		//if a solution with better cut_size was found AND all positive gains have been exhausted
		if(positive_gains_exhausted && cut_size <= min_cut_size)
		{
			saveBestSolution();
		}
	} 	
		    	
	//revert back to the best solution that was found this pass
	recallBestSolution();
	
	if(starting_cut_size == cut_size) //if no improvement was made i.e. stuck in local minimum
		return false;
	else return true;
}

//select the next base_cell
//if moving the cell would imbalance the partitioning, pick another cell
//finally, move the base_cell to the other partition
bool performNextMove()
{	
	if(gain_bucket[0].size() == 0 && gain_bucket[1].size() == 0) //if no more cell in gain bucket, end the pass
		return false;
		
	CELL* base_cell;	
	multimap <int, CELL*>::iterator itr [2]; //iterator for the gain buckets	
	
	if(gain_bucket[0].size() != 0)
		itr[0] = --gain_bucket[0].end(); //pointer to highest gain cell in partition 0
	if(gain_bucket[1].size() != 0)
		itr[1] = --gain_bucket[1].end();
		
	int from_partition;
			
	//if there is nothing left in one of the gain_buckets, use the other bucket
	if(gain_bucket[0].size() == 0) //if no cells in bucket 0
		from_partition = 1;
	else if(gain_bucket[1].size() == 0) //if no cells in bucket 1
		from_partition = 0;
	else from_partition = base_partition;
	
	base_partition = !base_partition; //pick from alternating partitions to make large nets dead faster		
	
	base_cell = itr[from_partition]->second; //pick the highest gain cell as default 

//if the highest gain cell would cause imbalance, pick another cell that doesn't imbalance
if(!cellCanMove(base_cell))
{
	//pick the highest gain cell from the other partition
	from_partition = !from_partition;
	base_cell = itr[from_partition]->second;
		
	if(!cellCanMove(base_cell)) //if the highest gain cell in BOTH partitions would cause an imbalance, no legal moves. End the pass
	{
		printf("No legal moves remain. Both would cause imbalance!\n");
		return false;
	}
}

	(*base_cell).changePartition();
	
	//return true indicating a successful move
	return true;
}

//check if the indicated cell can legally be moved to the other partition without causing imbalance
bool cellCanMove(CELL* base_cell)
{
	if(base_cell == NULL) { printf("base_cell is NULL!\n"); return false;}
	bool floor = (float)(partition_area[base_cell->partition] - base_cell->area) / (float)(partition_area[0] + partition_area[1]) >= floor_ratio;
	return floor;
}

void printResult(){
	cout << "******************************************" << endl;
	cout << "The result is:" << endl;
	for(size_t i = 0; i < cells.size(); i++){
		cout << "id:" << cells[i]->cell_id << " group:" << cells[i]->partition << " gain:" << cells[i]->gain << endl;
	}
	cout << "The final ";
	findCutsize();
}	

int main(){
	readInput();
	initialPartition();
	
	min_cut_size = cut_size; //initialize min_cut_size
	
	computeGains();	
	FMPartitionPass();
	printResult();
	return 0;
}
