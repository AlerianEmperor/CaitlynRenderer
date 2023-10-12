#ifndef _SORT_H_
#define _SORT_H_
#include "FlatNode.h"

//std:sort cause stack overflow if array is big, so I have to write a quicksort version
void swap(Reference* a, Reference* b)
{
	Reference t = *a;
	*a = *b;
	*b = t;
}

int partition(vector<Reference>& arr, int l, int h, int dim)
{
	Reference x = arr[h];
	int i = (l - 1);

	//float ca = r1.box.c()[dim];
	//float cb = r2.box.c()[dim];

	//return ca < cb;
	//return (ca < cb) || (ca == cb && r1.id < r2.id);

	for (int j = l; j <= h - 1; j++) 
	{
		float ca = arr[i].box.c()[dim];
		float cb = x.box.c()[dim];
		if(ca < cb || (abs(ca - cb) <= 0.0001f && arr[i].id < x.id))
		{
			i++;
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[h]);
	return (i + 1);
}

void quick_sort_iterative(vector<Reference>& arr, int l, int h, int dim)
{
	// Create an auxiliary stack
	vector<int> stack(h - l + 1);

	// initialize top of stack
	int top = -1;

	// push initial values of l and h to stack
	stack[++top] = l;
	stack[++top] = h;

	// Keep popping from stack while is not empty
	while (top >= 0) {
		// Pop h and l
		h = stack[top--];
		l = stack[top--];

		// Set pivot element at its correct position
		// in sorted array
		int p = partition(arr, l, h, dim);

		// If there are elements on left side of pivot,
		// then push left side to stack
		if (p - 1 > l) {
			stack[++top] = l;
			stack[++top] = p - 1;
		}

		// If there are elements on right side of pivot,
		// then push right side to stack
		if (p + 1 < h) {
			stack[++top] = p + 1;
			stack[++top] = h;
		}
	}
	vector<int>().swap(stack);
}

#endif // !_SORT_H_

