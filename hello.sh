#!/bin/bash          

#echo "Hello World"

#HELLO is the variable that holds the value 'Hello Vijay'
HELLO="Hello Vijay"

#echo is the command to print out the value of the variable HELLO
echo $HELLO

#command substitution
myfiles=$( ls ~/ | wc -l )
echo $myfiles

#variable a,b
a=11
b=11
echo $a, $b

#arithmetic operators
c=$(( a+b ))
echo "a+b" is $c
echo "a-b" is $(( a-b ))
echo "a*b" is $(( a*b ))

<<commentall
#user input, get two numbers
echo "Enter your first name : "
read firstname
echo "Enter your last name : "
read lastname
echo "Your Full Name is : $firstname $lastname"
commentall

#conditionals
if [ $b -gt $a ]
	then 
		echo "$b is the larger number"
elif [ $b -lt $a ]
	then
		echo "$a is the larger number"
	else
		echo "$a and $b are equal"
fi

#while loop
number=1
while [ $number -le 4 ]
	do
	echo $number
	((number++))
	done

#for loop list
aligners='bwa star bowtie tophat trinity'
for item in $aligners
	do
	echo $item
	done

#for loop range
for number in {1..10..1}
	do
	echo $number
	done

#a function that converts celcius to fahrenheit
function ctof 
	{
	ft=$(( $1 * (9/5) +  32 ))
	echo $ft
	}

#call the function
ctof 37



#the variable dna holds the value of a oligonucleotide
dna="atgcagtacagatataagcacac"


