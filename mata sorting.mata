mata
function partition(arr, left, right) {
	i = left
	j = right
	pivot = arr[(left + right) / 2]
	while (i <= j) {
		while (arr[i] < pivot) i++
        while (arr[j] > pivot) j--
		if (i <= j) {
			tmp = arr[i]
			arr[i] = arr[j]
			arr[j] = tmp
			i++
			j--
		}
    }
	return(i)
}

function quicksort(arr, left, right) {
	index = partition(arr, left, right)
	if (left < index - 1) quicksort(arr, left, index - 1)
	if (index < right) quicksort(arr, index, right)
}

function quickie(real colvector x){

	if (rows(x)<2){
		return(x)
	}
	else{
		el = x[1]
		part=x[2..rows(x)]
		v1 = select(part, part[.,1]:<=el)
		v2 = select(part, part[.,1]:>el)
		
		if (rows(v1)>0) v1 = quickie(v1)
		if (rows(v2)>0)	v2 = quickie(v2)
		if (rows(v1)>0 & rows(v2)>0) return(v1\el\v2)
		if (rows(v1)>0 & rows(v2)==0) return(v1\el)
		if (rows(v1)==0 & rows(v2)>0) return(el\v2)
	}
}

timer_clear()
x=rnormal(100000,1,10,5)

timer_on(2)
p=sort(x,1)
timer_off(2)

timer_on(1)
o=quickie(x)
timer_off(1)

timer_on(3)
quicksort(x,1,rows(x))
timer_off(3)
timer()
