// A basic implementation of insertion sort. It sorts an array a with n
// integers.
void isort(int *a, const int n)
{
	for (int i = 1; i < n; i++) {
		for (int j = i; j > 0 && a[j] < a[j-1]; j--) {
			int tmp = a[j];
			a[j] = a[j-1];
			a[j-1] = tmp;
		}
	}
}

