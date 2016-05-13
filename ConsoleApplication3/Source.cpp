#include <iostream>
#include <string>
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <list>


struct scoringMatrix
{
	char rowChar[24];
	char colChar[24];
	int **score;
};

class NeedlemanWunsch
{
private:
	int flagCheck = 0; //FLAG 1:CROSS,FLAG 2: LEFT,FLAG 3:UP
	int gap = 0;
	int match = 0;
	int unmatch = 0;
public:
	int matchCounter = 0;
	int positiveUnmatchCounter = 0;
	void nw(std::vector<char> fasta1, std::vector<char> fasta2, scoringMatrix matrix);
	int maxValue(int i, int j, int x);
};

int NeedlemanWunsch::maxValue(int i, int j, int x)
{
	if (i >= j && i >= x)
	{
		flagCheck = 1; //CROSS
		return i;
	}
	else if (j > x)
	{
		flagCheck = 2; // LEFT
		return j;
	}
	else
	{
		flagCheck = 3; // UP
		return x;
	}
}

void NeedlemanWunsch::nw(std::vector<char> fasta1, std::vector<char> fasta2, scoringMatrix matrix)
{
	std::vector<std::vector<int> > diagram;
	std::vector<std::vector<int> > trackDiagram;
	for (int i = 0; i <= fasta2.size()+1; i++)
	{
		std::vector<int> rows;
		for (int j = 0; j <= fasta1.size()+1; j++)
		{
			rows.push_back(j*i);
		}
		diagram.push_back(rows);
	}

	for (int i = 0; i <= fasta2.size() + 1; i++)
	{
		std::vector<int> rows;
		for (int j = 0; j <= fasta1.size() + 1; j++)
		{
			rows.push_back(j*i);
		}
		trackDiagram.push_back(rows);
	}
	gap = matrix.score[0][23]; //GAP SCORE
	for (int i = 0; i <= fasta1.size(); i++)
	{
		diagram[0][i] = 0;
	}
	for (int i = 0; i <= fasta2.size(); i++)
	{
		diagram[i][0] = 0;
	}
	int i = 0;
	int k = 0; //GAP counter
	int x = 0; //Similarity counter
	int j = 0;
	int valueCross = 0;
	int valueLeft = 0;
	int valueUp = 0;
	int temp = 0;
	for (i = 1; i <= fasta1.size(); i++)
	{
		for (j = 1; j <= fasta2.size(); j++)
		{
			int fasta1Value = 0;
			switch (fasta1[i - 1])
			{
			case 'A':fasta1Value = 0; break;
			case 'R':fasta1Value = 1; break;
			case 'N':fasta1Value = 2; break;
			case 'D':fasta1Value = 3; break;
			case 'C':fasta1Value = 4; break;
			case 'Q':fasta1Value = 5; break;
			case 'E':fasta1Value = 6; break;
			case 'G':fasta1Value = 7; break;
			case 'H':fasta1Value = 8; break;
			case 'I':fasta1Value = 9; break;
			case 'L':fasta1Value = 10; break;
			case 'K':fasta1Value = 11; break;
			case 'M':fasta1Value = 12; break;
			case 'F':fasta1Value = 13; break;
			case 'P':fasta1Value = 14; break;
			case 'S':fasta1Value = 15; break;
			case 'T':fasta1Value = 16; break;
			case 'W':fasta1Value = 17; break;
			case 'Y':fasta1Value = 18; break;
			case 'V':fasta1Value = 19; break;
			case 'B':fasta1Value = 20; break;
			case 'Z':fasta1Value = 21; break;
			case 'X':fasta1Value = 22; break;
			}
			int fasta2Value = 0;
			switch (fasta2[j - 1])
			{
			case 'A':fasta2Value = 0; break;
			case 'R':fasta2Value = 1; break;
			case 'N':fasta2Value = 2; break;
			case 'D':fasta2Value = 3; break;
			case 'C':fasta2Value = 4; break;
			case 'Q':fasta2Value = 5; break;
			case 'E':fasta2Value = 6; break;
			case 'G':fasta2Value = 7; break;
			case 'H':fasta2Value = 8; break;
			case 'I':fasta2Value = 9; break;
			case 'L':fasta2Value = 10; break;
			case 'K':fasta2Value = 11; break;
			case 'M':fasta2Value = 12; break;
			case 'F':fasta2Value = 13; break;
			case 'P':fasta2Value = 14; break;
			case 'S':fasta2Value = 15; break;
			case 'T':fasta2Value = 16; break;
			case 'W':fasta2Value = 17; break;
			case 'Y':fasta2Value = 18; break;
			case 'V':fasta2Value = 19; break;
			case 'B':fasta2Value = 20; break;
			case 'Z':fasta2Value = 21; break;
			case 'X':fasta2Value = 22; break;
			}
			valueCross = 0;
			valueLeft = 0;
			valueUp = 0;
			if (fasta1Value == fasta2Value)
			{
				match = matrix.score[fasta1Value][fasta2Value];
				valueCross = diagram[j - 1][i - 1] + match;
			}
			else
			{
				unmatch = matrix.score[fasta1Value][fasta2Value];
				valueCross = diagram[j - 1][i - 1] + unmatch;
			}
			valueLeft = diagram[j][i - 1] + gap;
			valueUp = diagram[j - 1][i] + gap;
			diagram[j][i] = maxValue(valueCross, valueLeft, valueUp);
			trackDiagram[j][i] = flagCheck;
		}
	}
	//TO SET INDEX OF ARRAY STARTUP
	std::cout << "Alignment score is: " << diagram[j-1][i-1] << '\n';
	i--; j--;
	std::list <char> alignedFasta1;
	std::list <char> alignedFasta2;
	std::list <char> alignment;
	int getNextFlag = flagCheck;
	while (i>0 && j>0) // RUN IN REVERSE ORDER THEN FLIP
	{
		switch (getNextFlag)
		{
		case 1:
			//MEANS WE HAVE MATCH OR UNMATCH
			alignedFasta1.push_front(fasta1[i-1]);
			alignedFasta2.push_front(fasta2[j-1]);
			if (alignedFasta1.front() == alignedFasta2.front()) //MATCH
			{
				matchCounter++;
				alignment.push_front('|');
			} 		
			else if (diagram[j][i] - diagram[j - 1][i - 1] < 0) //POSITIVE UNMATCH POINT COUNTER
			{
				positiveUnmatchCounter++;
				alignment.push_front(':');
			}
			else
			{
				alignment.push_front('.');
			}
			getNextFlag = trackDiagram[j - 1][i - 1];
			i--; j--;
			break;
		case 2:
			//MEANS WE NEED GAP AT 2ND FASTA
			alignedFasta1.push_front(fasta1[i-1]);
			alignedFasta2.push_front('-');
			getNextFlag = trackDiagram[j][i-1];
			alignment.push_front('-');
			i--;
			k++;
			break;
		case 3:
			//MEANS WE NEED GAP AT 1ST FASTA
			alignedFasta1.push_front('-');
			alignedFasta2.push_front(fasta2[j-1]);
			getNextFlag = trackDiagram[j - 1][i];
			alignment.push_front('-');
			j--;
			k++;
			break;
		}
		//PRINT OUT REMANING FASTA IF ONE REACHED TO ZERO
		if (i == 0)
		{
			for (; j == 1; j--)
			{
				alignedFasta1.push_front('-');
				alignedFasta2.push_front(fasta2[j - 1]);
				alignment.push_front('-');
				k++;
			}
		}
		else
		{
			for (; i == 1; i--)
			{
				alignedFasta1.push_front(fasta1[i - 1]);
				alignedFasta2.push_front('-');
				alignment.push_front('-');
				k++;
			}
		}
	}
	std::cout << "Identity is: " << matchCounter << std::endl;
	std::cout << "Similarity is: " << positiveUnmatchCounter + matchCounter << std::endl;
	std::cout << "Final length of alignment is: " << alignedFasta1.size() << std::endl;
	std::cout << "GAP counter is: "<< k << std::endl;

	for (std::list<char>::iterator it = alignedFasta1.begin(); it != alignedFasta1.end(); ++it)
	{
		std::cout << ' ' << *it;
	}
	std::cout << '\n';

	for (std::list<char>::iterator it = alignment.begin(); it != alignment.end(); ++it)
	{
		std::cout << ' ' << *it;
	}
	std::cout << '\n';

	for (std::list<char>::iterator it = alignedFasta2.begin(); it != alignedFasta2.end(); ++it)
	{
		std::cout << ' ' << *it;
	}
	std::cout << '\n';
}

int main(int argc, char **argv)
{
	if (argc <= 2)
	{
		std::cout << "Need more arguments to run" << std::endl;
		return 0;
	}

	std::ifstream fasta1(argv[1]);
	std::ifstream fasta2(argv[2]);
	if (!fasta1.good() || !fasta2.good())
	{
		std::cout << "Error while opening fasta sequences" << std::endl;
	}
	std::string lastline, identifierFasta1, identifierFasta2;
	std::vector <char> fasta1Actual, fasta2Actual;
	while (std::getline(fasta1, lastline))
	{
		if (lastline[0] == '>')
		{
			identifierFasta1 = lastline;
		}
		else
			for (int i = 0; i < lastline.length(); i++)
				fasta1Actual.push_back(lastline[i]);
	}
	while (std::getline(fasta2, lastline))
	{
		if (lastline[0] == '>')
		{
			identifierFasta2 = lastline;
		}
		else
			for (int i = 0; i < lastline.length(); i++)
				fasta2Actual.push_back(lastline[i]);
	}
	bool firstActualLine = true;
	std::string blosumArray[32];
	int lineNo = 0;
	std::ifstream myFileName("BLOSUM62.txt");
	scoringMatrix getScores;
	getScores.score = new int*[24];
	for (int i = 0; i < 24; i++)
	{
		getScores.score[i] = new int[24];
		for (int j = 0; j < 24; j++)
			getScores.score[i][j] = 0;
	}
	while (std::getline(myFileName, blosumArray[lineNo]))
	{
		if (blosumArray[lineNo][0] == '#')
		{
			lineNo++;
		}
		else
		{
			if (firstActualLine)
			{
				blosumArray[lineNo].erase(std::remove_if(blosumArray[lineNo].begin(), blosumArray[lineNo].end(), ::isspace), blosumArray[lineNo].end());
				for (int i = 0; i < blosumArray[lineNo].length(); i++)
					getScores.rowChar[i] = blosumArray[lineNo].at(i);
				firstActualLine = false;
			}

			else
			{
				int counter = 0;
				for (int i = 0; i < blosumArray[lineNo].length(); i++)
				{
					if (i == 0)
					{
						getScores.colChar[lineNo - 7] = blosumArray[lineNo][i];
						counter = 0;
					}

					else
					{
						if (counter < 24)
						{
							char temp[3];
							for (int j = 0; j < 3; j++) //READ 3 CHARACTERS STEP BY STEP
							{
								temp[j] = blosumArray[lineNo][3 * i - 2 + j];
							}

							if (temp[0] == ' ' && temp[1] == ' ') //CHECK POSITIVE VALUES LESS THAN 10
							{
								getScores.score[counter][lineNo - 7] = (temp[2] - '0');
							}

							else if (temp[0] == ' '&& temp[1] == '-') //CHECK NEGATIVE VALUES LESS THAN 10
							{
								getScores.score[counter][lineNo - 7] = -1 * (temp[2] - '0');
							}

							else if (temp[0] == ' ' && temp[1] != '-' && temp[1] != ' ') // CHECK POSITIVE VALUES GREATER THAN 10
							{
								getScores.score[counter][lineNo - 7] = (10 * (temp[1] - '0')) + (temp[2] - '0');
							}

							else if (temp[0] == '-' && temp[1] != ' ' && temp[1] != '-') // CHECK NEGATIVE VALUES GREATER THAN 10
							{
								getScores.score[counter][lineNo - 7] = (-10 * (temp[1] - '0')) + (temp[2] - '0');
							}
							counter++;
						}

					}
				}
			}
			lineNo++;
		}
	}
	NeedlemanWunsch calculateAlignment;
	calculateAlignment.nw(fasta1Actual, fasta2Actual, getScores);
	return 0;
}