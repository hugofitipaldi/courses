#Quiz: Scripting with Raw Input
#Quiz: Generate Messages
#Imagine you're a teacher who needs to send a message to each of your students reminding them of their missing assignments and grade in the class. You have each of their names, number of missing assignments, and grades on a spreadsheet and just have to insert them into placeholders in this message you came up with:

#Hi [insert student name],
#This is a reminder that you have [insert number of missing assignments] assignments left to submit before you can graduate. Your current grade is [insert current grade] and can increase to [insert potential grade] if you submit all assignments before the due date.

#You can just copy and paste this message to each student and manually insert the appropriate values each time, but instead you're going to write a program that does this for you. 
#Write a script that does the following:

#Ask for user input 3 times. Once for a list of names, once for a list of missing assignment counts, and once for a list of grades. Use this input to create lists for names, assignments, and grades.
#Use a loop to print the message for each student with the correct values. The potential grade is simply the current grade added to two times the number of missing assignments.

names =  input("Enter names separated by commas: ")# get and process input for a list of names
assignments = input("Enter assignments separated by commas: ").split(",")  # get and process input for a list of the number of assignments
grades = input("Enter grades separated by commas: ").split(",")  # get and process input for a list of grades

## message string to be used for each student
## HINT: use .format() with this string in your for loop
message = "Hi {},\n\nThis is a reminder that you have {} assignments left to \
submit before you can graduate. You're current grade is {} and can increase \
to {} if you submit all assignments before the due date.\n\n"

## write a for loop that iterates through each set of names, assignments, and grades to print each student's message

for name, assignment, grade in zip(names.split(","), assignments, grades):
    print(message.format(name, assignment, grade, int(grade) + int(assignment) * 2))