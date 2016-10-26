/* -------------------------------------------------------------------------- */

#include <iostream>
using namespace std;
#include <time.h>		// Used for clock().
#include <stdlib.h>		// Used for rand() and srand().

/* -------------------------------------------------------------------------- */

size_t	NumComparisons;
size_t	NumSwaps;

/* -------------------------------------------------------------------------- */

//
//	Minimise changes in the rest of the code, if we want to
//	change the basic data type (e.g. float or double).
//
typedef int TItem;

//
//	Swaps the values of the two arguments passed to it.
//
void Swap(TItem *Arg1, TItem *Arg2)
{
	TItem	Tmp = *Arg1;
	*Arg1 = *Arg2;
	*Arg2 = Tmp;
}	// Swap

//
//	Displays the contents of the array in the parameter.
//
void PrintArray(TItem Array[], size_t Size)
{
	size_t	i;


	cout << "Array contents:\n---------------\n";
	for ( i = 0;   i < ( Size - 1 );   i++ ) cout << "[" << i << "]=" << Array[i] << "  ";
	cout << "[" << i << "]=" << Array[i] << endl << endl;
}	// PrintArray

//
//	Displays the contents of a specified range of an array.
//
void PrintArray(TItem Array[], size_t Lo, size_t Hi)
{
	size_t	i;


	cout << "Array contents:\n---------------\n";
	for ( i = Lo;   i < Hi;   i++ ) cout << "[" << i << "]=" << Array[i] << "  ";
	cout << "[" << i << "]=" << Array[i] << endl << endl;
}	// PrintArray

//
//	Copies one array to another.
//
void CopyArray(TItem DestArray[], TItem SrcArray[], size_t Size)
{
	for ( size_t i = 0;   i < Size;   i++ ) DestArray[i] = SrcArray[i];
}	// CopyArray

//
//	Fills the array with random values, given
//	a seed for the random number generator.
//
void RandomiseArray(TItem Array[], size_t Size, unsigned int Seed)
{
	//	Make sure we always get the same random number sequence.
	srand(Seed);

	//	Fill the array with random values.
	for ( size_t i = 0;   i < Size;   i++ ) Array[i] = rand();
}	// RandomizeArray

//
//	Returns true if the array is sorted, else false.
//
bool IsArraySorted(TItem Array[], size_t Size)
{
	size_t	i = 1;
	bool	Sorted = true;


	while ( ( i < Size ) && Sorted )
	{
		Sorted = ( Array[i - 1] <= Array[i] );
		i++;
	}	// while

	if ( !Sorted ) cout << "IsArraySorted: i = " << i - 1 << endl;
/*
if ( !Sorted )
{
//	size_t	tmp;
//	for ( tmp = ( i - 1 );   tmp < ( i + 2 );   tmp++ ) cout << Array[tmp] << endl;
	cout << Array[i - 2] << endl;
	cout << Array[i - 1] << endl;
	cout << Array[i - 0] << endl;
}
*/
	return ( Sorted );
}	// IsArraySorted

/* -------------------------------------------------------------------------- */
//
//	BUBBLE SORT
//

void BubbleSort(TItem Array[], size_t Size)
{
	size_t	i, j;


	for ( i = 0;   i < ( Size - 1 );   i++ )
	{
		for ( j = ( Size - 1 );   j > i;   j-- )
		{
			NumComparisons++;
			if ( Array[j] < Array[j - 1] )
			{
				NumSwaps++;
				Swap(&(Array[j]), &(Array[j - 1]));
			}	// if
		}	// for
	}	// for
}	// BubbleSort

void BubbleSort(TItem Array[], size_t Lo, size_t Hi)
{
	size_t	i, j;


	for ( i = Lo;   i < ( Hi - 1 );   i++ )
	{
		for ( j = ( Hi - 1 );   j > i;   j-- )
		{
			if ( Array[j] < Array[j - 1] ) Swap(&(Array[j]), &(Array[j - 1]));
		}	// for
	}	// for
}	// BubbleSort

void BubbleSortOpt(TItem Array[], size_t Size)
{
	size_t	i, n = Size;
	bool	Swapped;


	//	Slightly optimised bubble sort.
	do
	{
		Swapped = false;

		for ( i = 0;   i < ( n - 1 );   i++ )
		{
			NumComparisons++;
			if ( Array[i] > Array[i + 1] )
			{
				NumSwaps++;
				Swap(&(Array[i]), &(Array[i + 1]));
				Swapped = true;
			}	// if
		}	// for

		n--;
	}
	while ( Swapped );
}	// BubbleSortOpt

void BubbleSortOpt(TItem Array[], size_t Lo, size_t Hi)
{
	size_t	i, n = ( Hi - Lo + 1 );
	bool	Swapped;


	//	Slightly optimised bubble sort.
	do
	{
		Swapped = false;

		for ( i = 0;   i < ( n - 1 );   i++ )
		{
			if ( Array[i] > Array[i + 1] )
			{
				Swap(&(Array[i]), &(Array[i + 1]));
				Swapped = true;
			}	// if
		}	// for

		n--;
	}
	while ( Swapped );
}	// BubbleSortOpt

/* -------------------------------------------------------------------------- */
//
//	INSERTION SORT
//

void InsertionSort(TItem Array[], long int Size)
{
	TItem		Key;
	long int	i, j;


	for ( i = 1;   i < Size;   i++ )
	{
		Key = Array[i];
		j = ( i - 1 );

		while ( ( j >= 0 ) && ( Key < Array[j] ) )
		{
			Array[j + 1] = Array[j];
			j--;
		}	// while

		Array[j + 1] = Key;
	}	// for
}	// InsertionSort

void InsertionSort(TItem Array[], long int Lo, long int Hi)
{
	TItem		Key;
	long int	i, j;


	for ( i = ( Lo + 1 );   i < ( Hi - Lo + 1 );   i++ )
	{
		Key = Array[i];
		j = ( i - 1 );

		while ( ( j >= 0 ) && ( Key < Array[j] ) )
		{
			Array[j + 1] = Array[j];
			j--;
		}	// while

		Array[j + 1] = Key;
	}	// for
}	// InsertionSort

/* -------------------------------------------------------------------------- */
//
//	SELECTION SORT
//

void SelectionSort(TItem Array[], size_t Size)
{
	size_t	i, j, Min;


	for ( j = 0;   j < Size;   j++ )
	{
		//	Find the minimum element of unsorted sub-array Array[j ... Size-1].
		Min = j;	//	Assume the minimum element is the first one.

		//	Find the smallest element after j.
		for ( i = ( j + 1 );   i < Size;   i++ )
		{
			//	If there is a smaller element, this is the new minimum.
			if ( Array[i] < Array[Min] ) Min = i;	// Found new min.
		}	// for

		//	Min is the index of the minimum element.
		//	Swap it with the current position.
		if ( Min != j ) Swap(&(Array[j]), &(Array[Min]));
	} // for
}	// SelectionSort

void DoubleSelectionSort(TItem Array[], size_t Size)
{
	size_t	MinIndex, MaxIndex, L = 0, R = ( Size - 1 ), i;
	bool	Swapped;


	if ( Size < 2 ) return;		// Nothing to sort.

	do
	{
		MinIndex = MaxIndex = L;	// Start of sub-array.
		Swapped = false;			// Haven't swapped any elements yet.

		//	Find the positions of the minimum and maximum elements.
		for ( i = ( L + 1 );   i <= R;   i++ )
		{
			if ( Array[i] < Array[MinIndex] ) MinIndex = i;		// New min index.
			if ( Array[i] > Array[MaxIndex] ) MaxIndex = i;		// New max index.
		}	// for

		// Place minimum value in the first position of the current sub-array.
		if ( L != MinIndex )
		{
			Swap(&(Array[L]), &(Array[MinIndex]));
			Swapped = true;
		}	// if

		//	Check if we swapped L with the MaxIndex in the previous IF statement.
		if ( L == MaxIndex ) MaxIndex = MinIndex;

		// Place maximum value in the last position of the current sub-array.
		if ( R != MaxIndex )
		{
			Swap(&(Array[R]), &(Array[MaxIndex]));
			Swapped = true;
		}	// if

		L++;	// Move the left sub-index one element to the right.
		R--;	// Move the right sub-index one element to the left.
	}
	while ( Swapped && ( L < R ) );
}	// DoubleSelectionSort

void DoubleSelectionSort(TItem Array[], size_t Lo, size_t Hi)
{
	size_t	MinIndex, MaxIndex, L = Lo, R = Hi, i;


	while ( L < R )
	{
		MinIndex = MaxIndex = L;	// Start of sub-array.

		//	Find the positions of the minimum and maximum elements.
		for ( i = ( L + 1 );   i <= R;   i++ )
		{
			if ( Array[i] < Array[MinIndex] ) MinIndex = i;		// New min index.
			if ( Array[i] > Array[MaxIndex] ) MaxIndex = i;		// New max index.
		}	// for

		// Place minimum value in the first position of the current sub-array.
		if ( L != MinIndex ) Swap(&(Array[L]), &(Array[MinIndex]));

		//	Check if we swapped L with the MaxIndex in the previous IF statement.
		if ( L == MaxIndex ) MaxIndex = MinIndex;

		// Place maximum value in the last position of the current sub-array.
		if ( R != MaxIndex ) Swap(&(Array[R]), &(Array[MaxIndex]));

		L++;	// Move the left sub-index one element to the right.
		R--;	// Move the right sub-index one element to the left.
	}	// while
}	// DoubleSelectionSort

/* -------------------------------------------------------------------------- */
//
//	QUICK SORT
//

void QuickSort(TItem Array[], size_t L, size_t R)
{
	TItem	Key;
	size_t	i, j, k;	// k = Pivot index.



	if ( L < R )	// Check if already sorted. --- Non-degenerate case.
	{
		//k = ( ( L + R ) / 2 );	// Select initial pivot position.
		k = L + ( ( R - L ) / 2 );	// see Wikipedia on quick sort.
		Swap(&(Array[L]), &(Array[k]));
		Key = Array[L];
		i = ( L + 1 );
		j = R;

		while ( i <= j )
		{
			while ( ( i <= R ) && ( Array[i] <= Key ) ) i++;
			while ( ( j >= L ) && ( Array[j] > Key ) ) j--;
			if( i < j ) Swap(&(Array[i]), &(Array[j]));
		}	// while

		//	Swap two elements.
		Swap(&(Array[L]), &(Array[j]));

		//	Recursively sort the two sub-lists.
		QuickSort(Array, L, j - 1);
		QuickSort(Array, j + 1, R);
	}	// if
}	// QuickSort
/*
void QuickSortHybrid(TItem Array[], int L, int R)
{
	TItem	Key;
	int		i, j, k;	// k = Pivot index.


	if ( L >= R ) return;	// Nothing to sort.

	if ( ( R - L ) <= 5 )
	{
		//	Degenerate case: array partition too small for recursion.
		InsertionSort(Array, L, R);
	}
	else
	{
		//	Non-degenerate case.
		k = L + ( ( R - L ) / 2 );	// Select initial pivot position without overflow.
		Swap(&(Array[L]), &(Array[k]));
		Key = Array[L];
		i = ( L + 1 );
		j = R;

		while ( i <= j )
		{
			while ( ( i <= R ) && ( Array[i] <= Key ) ) i++;
			while ( ( j >= L ) && ( Array[j] > Key ) ) j--;
			if( i < j ) Swap(&(Array[i]), &(Array[j]));
		}	// while

		//	Swap two elements.
		Swap(&(Array[L]), &(Array[j]));

		//	Recursively sort the two sub-lists.
		QuickSortHybrid(Array, L, j - 1);
		QuickSortHybrid(Array, j + 1, R);
	}	// if
}	// QuickSortHybrid
*/
/*
http://www.algolist.net/Algorithms/Sorting/Quicksort
http://www.cs.wcupa.edu/~rkline/DS/fast-sorts.html
*/
/*
void QuickSortHybrid2(TItem Array[], size_t Lo, size_t Hi)
{
	TItem	Pivot = Array[Lo + ( ( Hi - Lo ) / 2 )];
	size_t	i = Lo, j = Hi;


	//	Partition.
	while ( i <= j )
	{
		while ( Array[i] < Pivot ) i++;
		while ( Array[j] > Pivot ) j--;
		if ( i <= j ) Swap(&(Array[i++]), &(Array[j--]));
	}	// while

	//	Recursion.
	if ( Lo < j )
	{
//		PrintArray(Array, Lo, j);
//		cout << "Sort: Lo = " << Lo << ", j = " << j << endl;
		QuickSortHybrid2(Array, Lo, j);
	}
	if ( i < Hi )
	{
//		PrintArray(Array, i, Hi);
//		cout << "Sort: i = " << i << ", Hi = " << Hi << endl;
		QuickSortHybrid2(Array, i, Hi);
	}
}	// QuickSortHybrid2
*/
/*
void QSort(TItem arr[], int left, int right)
{
	int		i = left, j = right;
	int		tmp;
	int		pivot = arr[(left + right) / 2];


	// partition
	while (i <= j)
	{
		while (arr[i] < pivot) i++;
		while (arr[j] > pivot) j--;
		if (i <= j)
		{
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;
			j--;
		}
	}

	// recursion
	if (left < j) QSort(arr, left, j);
	if (i < right) QSort(arr, i, right);
}	// QSort
*/

void QSort(TItem Array[], size_t left, size_t right)
{
	size_t	i = left, j = right;
	TItem	pivot = Array[left + ( right - left ) / 2];


	//	1. Partition
	while ( i <= j )
	{
		while (Array[i] < pivot) i++;
		while (Array[j] > pivot) j--;
		if ( i <= j )
		{
			Swap(&(Array[i]), &(Array[j]));
			i++;
			j--;
		}	// if
	}	// while

	//	2. Recursion
	if ( left < j ) QSort(Array, left, j);
	if ( i < right ) QSort(Array, i, right);
}	// QSort

/* -------------------------------------------------------------------------- */
//
//	MERGE SORT
//

void MergeSort(TItem A[], int Temp[], int L, int R)
{
	//	Condition that terminates recursion.  The degenerate case.
	if ( L == R ) return;

	//	Midpoint selection.
	//int Mid = ( ( L + R ) / 2 );
	int Mid = L + ( ( R - L ) / 2 );	// see Wikipedia on quick sort.

	//	Recursively split array in half.
	MergeSort(A, Temp, L, Mid);
	MergeSort(A, Temp, Mid + 1, R);

	//	Start the process of merging.
	int		k = L, i = L, j = Mid + 1;

	while ( ( i <= Mid ) && ( j <= R ) )
	{
		if ( A[i] < A[j] )
		{
			Temp[k] = A[i];
			i++;
		}
		else
		{
			Temp[k] = A[j];
			j++;
		}	// if
		k++;
	}	// while

	//	Copy remaining elements at the end of Temp.
	while ( i <= Mid )
	{
		Temp[k] = A[i];
		k++;
		i++;
	}	// while

	while ( j <= R )
	{
		Temp[k] = A[j];
		k++;
		j++;
	}	// while

	//	Copy everything back to A.
	for ( i = L;   i <= R;   i++ ) A[i] = Temp[i];
}	// MergeSort

/* -------------------------------------------------------------------------- */

#define		MAX		100000
//#define		MAX		10000000
//#define		MAX		20

int main()
{
	TItem		*Arr, *Tmp;
	clock_t		Start, Stop;
	double		Sec;


	//	Attempt to allocate memory for array.
	Arr = new TItem [MAX];
	if ( Arr == NULL )
	{
		cout << "Failed to allocate main array.  Aborting" << endl;
		exit(1);
	}	// if

	Tmp = new TItem[MAX];
	if ( Tmp == NULL )
	{
		delete [] Arr;
		cout << "Failed to allocate auxiliary array.  Aborting" << endl;
		exit(2);
	}	// if

	//	Test BubbleSort.
	NumComparisons = 0;
	NumSwaps = 0;
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
	cout << "Sorting..." << endl;
	Start = clock();
	BubbleSort(Arr, MAX);
	Stop = clock();
	Sec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC );
	cout << "Time to execute Bubble Sort: " << Sec << " sec." << endl << endl;
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl;// << endl;
	cout << "Comparisons = " << NumComparisons << endl;
	cout << "      Swaps = " << NumSwaps << endl << endl;

	/*//	Test BubbleSortOpt.
	NumComparisons = 0;
	NumSwaps = 0;
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
	cout << "Sorting..." << endl;
	Start = clock();
	BubbleSortOpt(Arr, MAX);
	Stop = clock();
	Sec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC );
	cout << "Time to execute Bubble Sort (Opt): " << Sec << " sec." << endl << endl;
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl;// << endl;
	cout << "Comparisons = " << NumComparisons << endl;
	cout << "      Swaps = " << NumSwaps << endl << endl;

	//	Test InsertionSort.
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
	cout << "Sorting..." << endl;
	Start = clock();
	InsertionSort(Arr, MAX);
	Stop = clock();
	Sec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC );
	cout << "Time to execute Insertion Sort: " << Sec << " sec." << endl;
//	PrintArray(Arr, MAX);
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl << endl;

	//	Test SelectionSort.
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
	cout << "Sorting..." << endl;
	Start = clock();
	SelectionSort(Arr, MAX);
	Stop = clock();
	Sec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC );
	cout << "Time to execute Selection Sort: " << Sec << " sec." << endl << endl;
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl << endl;

	//	Test DoubleSelectionSort.
//Itr = 0;
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
//PrintArray(Arr, MAX);
	cout << "Sorting..." << endl;
	Start = clock();
//	DoubleSelectionSort(Arr, MAX);
	DoubleSelectionSort(Arr, 0, MAX - 1);
	Stop = clock();
	Sec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC );
//PrintArray(Arr, MAX);
	cout << "Time to execute Double Selection Sort: " << Sec << " sec." << endl << endl;
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl << endl;
//	cout << "ITR = " << Itr << endl;


	//	Test QuickSort.
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
	cout << "Sorting..." << endl;
	Start = clock();
	QuickSort(Arr, 0, MAX - 1);
	Stop = clock();
	Sec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC );
	cout << "Time to execute Quick Sort: " << Sec << " sec." << endl;
//	PrintArray(Arr, MAX);
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl << endl;
/*
	//	Test QuickSortOpt.
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
	cout << "Sorting..." << endl;
	Start = clock();
	QuickSortHybrid(Arr, 0, MAX - 1);
	Stop = clock();
	Msec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC * 1000.0 );
	cout << "Time to execute Quick Sort (Hybrid): " << Msec << " msec." << endl;
//	PrintArray(Arr, MAX);
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl << endl;
*/
/*
	//	Test QuickSortHybrid2.
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
	cout << "Sorting..." << endl;
	Start = clock();
	QuickSortHybrid2(Arr, 0, MAX - 1);
	Stop = clock();
	Msec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC * 1000.0 );
	cout << "Time to execute Quick Sort (Hybrid2): " << Msec << " msec." << endl;
//	PrintArray(Arr, MAX);
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl << endl;
*/

	//	Test QSort.
	/*
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
	cout << "Sorting..." << endl;
	Start = clock();
	QSort(Arr, 0, MAX - 1);
	Stop = clock();
	Sec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC );
	cout << "Time to execute QSort: " << Sec << " sec." << endl;
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl << endl;

	//	Test MergeSort.
	cout << "Randomising array..." << endl;
	RandomiseArray(Arr, MAX, 1000);
	cout << "Sorting..." << endl;
	Start = clock();
	MergeSort(Arr, Tmp, 0, MAX - 1);
	Stop = clock();
	Sec = ( (double) ( Stop - Start ) / CLOCKS_PER_SEC );
	cout << "Time to execute Merge Sort: " << Sec << " sec." << endl;
//	PrintArray(Arr, MAX);
	cout << "IsArraySorted = " << IsArraySorted(Arr, MAX) << endl << endl; */

	//	Deallocate memory.
	delete [] Arr;
	delete [] Tmp;
	return 0;
}	// main

/* -------------------------------------------------------------------------- */
