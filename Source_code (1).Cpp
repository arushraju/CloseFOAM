/*Welcome to this software where we will be making a FLUID SIMULATOR, and this will be
done using finite difference method.
There are many parts of this project, and this is one of them. Actually in the first part of the project
essentially, I will be building this software which will take in a matrix of size 
lets say 1000*1000 pixel, and the matrix must always contain element from the set 0,1,2,3 only.
*/

/*When thinking about the libraries there are many options I have*/ 
//Standard libraries
#include<stdio.h> //Surely for input and output.
#include<stdlib.h> //This is for the memory allocation.
#include<memory.h> //This is for the space complexity
#include<time.h> //This is to keep track of the time complexity
#include<math.h> //This is for usign mathematical functions
#include<string.h>//This for the terminal commands
//For High Performance Computing (openmp)
#include<omp.h> //This is for parallel computing
//For Animations and terminal front end
#include<conio.h> // For options
#include<windows.h>//This is for terminal animation


//This is for the CSV file ouptut
FILE* stream_vorticity_velocity_fr; //This will store velocity at all the time step


bool file_output = true;//Do you want the file to be outputed in the form of csv?
bool plot = true; //This will automatically run the terminal plot command for python script

//This is for debugging
bool debugger = false;
bool strong_debugger = false;//This is for even more detailed debugging
bool computation_time = true; //This will show the computation time
bool terminate = false;//This will be used to terminate the program where ever required


/*This is the structure that will contain two 1D array(pointer).The first one will contain boundary matrix that
will be the matrux that will store the boundary condition, like the type of boundary, and the value matrx will
store the values of unknow quantity, like the colors intially*/
typedef struct input_data_1{
    int* boundary_matrix;
    double* value_matrix;
}input_1;

/*The same idea could be extended to 2D. Here the boundary conditions are as follows. The first matrix will contain
boundary matrix that will be 2D array. And finally the value matrix will be matrix with 2D array that will contain
all the values at where put into inout initially */
typedef struct input_data_2{
    int** boundary_matrix;
    double** value_matrix;
}input_2;


//This will be used to set the type of problem (Eg: Heat Equation)
int dimension;//Dimension of problem (1,2)
int problem_type;//This could go from (1 to 20)
int scheme_type;//Kind of scheme to be used (implicit(2) or explicit(1))
int scheme;      //Kind of scheme to be used, and this number would depend on the above two numbers.

// These are the dicretization points of the grid in different dimension in space and time
int grid_points_x;
int grid_points_y;
int time_steps;

// These will be the dimension of the doamin and time
double L_x;
double L_y;
double time_max;

/*This function will be used to solve the problem it should be called and it will make the move*/
void print_output_1(input_1 input)
{
    //Printing the output
    printf("\n");
    for(int i=0;i<grid_points_x;i++){printf("%d ",input.value_matrix[i]);};
    printf("\n");
    return;
}
void print_output_2(input_2 input)
{
    // Printing the values of output
    printf("\n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            printf("%lf\t",input.value_matrix[i][j]);
        }; printf("\n");};
    printf("\n");
    return;
}

/****************************************    HEAT EQUATION   ************************************ */

//Forward in time central in space scheme for one dimension (Explicit)
//So the first thing that needs to be done, is to make an array that will encode all the values of colors to numbers as intial condition
/*So here is the game plan. I surely know that the initial condition in this case has to be specified
and and it could be done by coloring the domain from red to blue in different shades and this could be encoded 
as numbers. So the plan has now changed. While solving this problem of heat eqaution, I will not ask for three
colors but the range of color, as if the contour is first colored and then the picture is aploaded.
As if the colors will then come to life with time*/
void FTCS_Heat_1(input_1 input,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3,double alpha)
{    
    //if debugging is required
    if(debugger){
        printf("Into the FTCS function\n");
        printf("Values stored inside the boundary matrix are: \n");
        for(int i=0;i<grid_points_x;i++){printf("%d ",input.boundary_matrix[i]);};
        printf("Values stored inside the values matrix are: \n");
        for(int i=0;i<grid_points_x;i++){printf("%lf ",input.value_matrix[i]);};
        printf("\n");
    }
    
    // This will be used for calculating the important numbers.
    double delta_x = L_x/grid_points_x;
    double delta_t = time_max/time_steps;
    
    //Clock starts
    clock_t start = clock();

    for(int time=0;time<time_steps;time++){
        if(debugger){printf("Time = %d\n",time);} 
        for(int space=0;space<grid_points_x;space++)
        {
            if(debugger)
            {
                printf("space = %d\n",space);
                printf("input.value_matrix[%d] = %lf\n",space,input.value_matrix[space]);
                printf("input.boundary_matrix[%d] = %d\n",space,input.boundary_matrix[space]);
            }

            //If the point is meant to be computed and updated
            //If the input is meant to be enforced by dirichlet_condition
            if(input.boundary_matrix[space]!=0){
                switch(input.boundary_matrix[space]) 
                {
                    case 1:
                        input.value_matrix[space] = dirichlet_condition_1;
                        break;
                    case 2:
                        input.value_matrix[space] = dirichlet_condition_2;
                        break;
                    case 3:
                        input.value_matrix[space] = dirichlet_condition_3;
                        break;
                    default:
                        printf("Unexpected entry in the boundary matrix\n");
                 return;
                }
            }

            else if(input.boundary_matrix[space]==0)
            {   
                if(debugger)
                {
                    printf("alpha = %lf|delta_t = %lf|(deltax^2) = %lf|r = %lf\n",alpha,delta_t,(pow(delta_x,2)),((alpha*delta_t)/(pow(delta_x,2))));
                    printf("input.space[%d]=%lf\n",space-1,input.value_matrix[space-1]);
                    printf("input.space[%d]=%lf\n",space+1,input.value_matrix[space+1]);
                    printf("input.space[%d]=%lf\n",space,input.value_matrix[space]);
                    printf("input.value_matrix[%d] + input.value_matrix[%d] - 2*input.value_matrix[%d]=%lf\n",space-1,space+1,space,input.value_matrix[space-1] + input.value_matrix[space+1] - 2*input.value_matrix[space]);   
                }
                
                input.value_matrix[space] = input.value_matrix[space] + ((alpha*delta_t)/(pow(delta_x,2)))*(input.value_matrix[space-1] + input.value_matrix[space+1] - 2*input.value_matrix[space]);
            }
            //Here we will print the output in the terminal
            if(debugger) printf("input.value_matrix[%d] = %lf\n",space,input.value_matrix[space]);
        }
        printf("\n");
    }
    
    //Clock ends
    clock_t end = clock();
    //Time taken
    if(computation_time){
        double cpu_time = (double)(end-start)/CLOCKS_PER_SEC;
        printf("Time taken : %f secs\n",cpu_time);
    }
    return;
};

//Dufort_Frankel_Scheme in one dimensional (Explict)
void Dufrot_Frankel_Heat_1(input_1 input,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3,double alpha){
    printf("Into the Dufrot Frankel Scheme function\n");
    double delta_x = L_x/grid_points_x;
    double delta_t = time_max/time_steps;
    double r = (alpha*delta_t)/(pow(delta_x,2));
    for(int time=0;time<time_steps;time++){
        printf("Time = %d\n",time);
        for(int space=0;space<grid_points_x;space++)
        {
            if(input.boundary_matrix[space]!=0){
                switch(input.boundary_matrix[space]) 
                {
                    case 1:
                        input.value_matrix[space] = dirichlet_condition_1;
                        break;
                    case 2:
                        input.value_matrix[space] = dirichlet_condition_2;
                        break;
                    case 3:
                        input.value_matrix[space] = dirichlet_condition_3;
                        break;
                    default:
                        printf("Unexpected entry in the boundary matrix\n");
                 return;
                }
            }
            
            else if(input.boundary_matrix[space]==0)
            {   
                //printf("alpha = %lf|delta_t = %lf|(deltax^2) = %lf|r = %lf\n",alpha,delta_t,(pow(delta_x,2)),((alpha*delta_t)/(pow(delta_x,2))));
                //printf("input.space[%d]=%lf\n",space-1,input.value_matrix[space-1]);
                //printf("input.space[%d]=%lf\n",space+1,input.value_matrix[space+1]);
                //printf("input.space[%d]=%lf\n",space,input.value_matrix[space]);
                //printf("input.value_matrix[%d] + input.value_matrix[%d] - 2*input.value_matrix[%d]=%lf\n",space-1,space+1,space,input.value_matrix[space-1] + input.value_matrix[space+1] - 2*input.value_matrix[space]);
                input.value_matrix[space] = ((1-2*r)/(1+2*r))*input.value_matrix[space] + (2/(1+2*r))*(input.value_matrix[space-1] + input.value_matrix[space+1]);
            }
            //Here we will print the output in the terminal
            printf("input.value_matrix[%d] = %lf\n",space,input.value_matrix[space]);
        }
        printf("\n");
    }
};

//Forward in time central in space scheme for two dimensional (Explicit)
void FTCS_Heat_2(input_2 input,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3,double alpha)
{
    printf("Into the FTCS function\n");
    double delta_x = L_x/grid_points_x;
    double delta_y = L_y/grid_points_y;
    double delta_t = time_max/time_steps;
    printf("time steps = %d\n",time_steps);

    double r_x = (alpha*delta_t)/(pow(delta_x,2));
    double r_y = (alpha*delta_t)/(pow(delta_y,2));

    for(int time=0;time<time_steps;time++){
        printf("Time = %d\n",time);
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0;i<grid_points_x;i++)
            {
                //If the input is meant to be enforced by dirichlet_condition
                if(input.boundary_matrix[i][j]!=0){
                    //printf("input.boundary_matrix[%d][%d] = %d\n",i,j,input.boundary_matrix[i][j]);
                    switch(input.boundary_matrix[i][j]) 
                    {
                        case 1:
                            input.value_matrix[i][j] = dirichlet_condition_1;
                            break;
                        case 2:
                            input.value_matrix[i][j] = dirichlet_condition_2;
                            break;
                        case 3:
                            input.value_matrix[i][j] = dirichlet_condition_3;
                            break;
                        default:
                            printf("Unexpected entry in the boundary matrix\n");
                    return;
                    }
                }

                else if(input.boundary_matrix[i][j]==0)
                {   
                    //printf("input.value_matrix[%d][%d] = %lf\n",i,j-1,input.value_matrix[i][j-1]);
                    input.value_matrix[i][j] = input.value_matrix[i][j] + 
                    (r_x)*(input.value_matrix[i+1][j] + input.value_matrix[i-1][j] - 2*input.value_matrix[i][j])+
                    (r_y)*(input.value_matrix[i][j+1] + input.value_matrix[i][j-1] - 2*input.value_matrix[i][j]);
                }
                //Here we will print the output in the terminal
                printf("%lf ",input.value_matrix[i][j]);
            }
            printf("\n");
        }
    }
}

//Dufort_Frankel_Scheme in two dimensional (Explict)
void Dufrot_Frankel_Heat_2(input_2 input,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3,double alpha)
{
    return;
}

input_1 Thomas_Algorithm_Heat_1(input_1 input,double r,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    bool debugger = false;
    printf("Into the Thomas_Algorithm_Heat_1 function\n");
    // Formation of matrix A
    /*So the thing is if grid points are of number n, then the matrix formed will be (n-2)*(n-2) since two
    points will be at the boundary, and they need not be computed*/
    //Now interesting thing is in the case of Crank Nicolsen method, the matrix formed will be tri-diagonal matrix
    int n = grid_points_x-2; //Dimension of A & b
    /*This is a bit tricky business to deal with implicit functions and it is necessary for us to stay careful.So for that reason I ma intalising the elements with zero when forming it,. Therefore I am using calloc function rather than malloc*/
    double** A = (double**)malloc(sizeof(double*)*(n));
    //its takes in two argument, the number of elements and the size of each element. Interesting!
    for(int i=0;i<n;i++){A[i] = (double*)calloc(n,sizeof(double));}
    if(debugger){
        printf("Just for debugging purpose printing the matrix A\n");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
    }
    
    /*So here is the convention. The row is first index of A, and column is second index of A. And the top right elements will be the A[0][0], and bottom left element will be A[n-2][n-2]*/
    for(int row=0;row<n;row++)
    {
        if(row>0) A[row][row-1] = -1;
        A[row][row] = ((2+2*r)/r);
        if(row<n-1) A[row][row+1] = -1;
    }
    if(debugger){
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A\n");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
    }

    //Creating b matrix
    /*This is possible to create at a specific time. Since the input that we have can be used to calculate b*/
    double* b = (double*)malloc(sizeof(double)*n);
    for(int row=0;row<n;row++){b[row] = input.value_matrix[row]+((2-2*r)/r)*input.value_matrix[row+1]+input.value_matrix[row+2];}
    if(debugger){
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after initialization\n");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");
    }
    /*And now the A and b matrix will be created. And now we could use the row replacement method to apply the
    dirichlet boundary condition in the interior points of the domain*/
    /*Lets say that I impose the value of u_row to be eqaul to u_imposed, then I must replace A[row][row] = 1 and b[row] = u_imposed, and
    then replace all the other values of matrix A with zero in that row.*/
    for(int l=1;l<grid_points_x-1;l++)
    {
        if(input.boundary_matrix[l])
        {
            for(int k=0;k<n;k++) {A[l-1][k] = 0;}
            A[l-1][l-1] = 1;
            switch (input.boundary_matrix[l])
            {
                case 1:
                    b[l-1] = dirichlet_condition_1;
                    break;
                case 2:
                    b[l-1] = dirichlet_condition_2;
                    break;
                case 3:
                    b[l-1] = dirichlet_condition_3;
                    break;
                default:
                    break;
            }
        }
    }
    if(debugger){
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after imposing boundary condition\n");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after imposing boundary condition\n");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");
    }
    
    //And now starts the Thomas Algorithm
    
    //Forward Elimination
    double temp;
    for(int row=0;row<n-1;row++)
    {
        //First we make the diagonal element equal to 1
        temp= A[row][row];
        for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
        b[row] = b[row]/temp;
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after making diagonal element equal to 1 at row = %d\n",row);
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after making diagonal element equal to 1 at row = %d\n",row);
            for(int i=0;i<n;i++){printf("%lf ",b[i]);}
            printf("\n");
        }
        // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
        temp = A[row+1][row];
        for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
        b[row+1] = b[row+1] - temp*b[row];
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after eliminating the element\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after eliminating the element\n");
            for(int i=0;i<n;i++){printf("%lf ",b[i]);}
            printf("\n");
        }
    }
    if(debugger){
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after forward elimination\n");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after forward elimination\n");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");
    }
    
    //Backward substitution
    input.value_matrix[n] = b[n-1]/A[n-1][n-1];
    for(int row=n-1;row>=1;row--)
    {
        input.value_matrix[row] = b[row-1]-A[row-1][row]*input.value_matrix[row+1];
    }
    if(debugger){
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after backward substitution\n");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after backward substitution\n");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the value matrix after backward substitution\n");
        for(int i=0;i<grid_points_x;i++){printf("%lf ",input.value_matrix[i]);}
        printf("\n");
    }
    //And this back substitution will finally complete all the problem, and we will be able to update all the inut values
    return input;
}

// Crank Nicolson Scheme in one dimensional (Implicit)
void Crank_Nicolson_scheme_Heat_1(input_1 input,int dirichlet_condition_1,int dirichlet_condition_2,int dirichlet_condition_3,double alpha)
{
    //Now this is a bit tricky. i will have to create the matrix A that will be present in Ax=b.
    //What I have with me is matrix that contain all the elements in value matrix and boundary matrix.
    printf("Into the Crank Nicolson_1 function\n");

    double delta_x = L_x/grid_points_x;
    double delta_t = time_max/time_steps;
    double r = (alpha*delta_t)/(pow(delta_x,2));
    for(int time=0;time<time_steps;time++){
        input = Thomas_Algorithm_Heat_1(input,r,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
        for(int x=0;x<grid_points_x;x++) printf("%lf\n",x,input.value_matrix[x]);
    }
    return;
}

/*This will be done later after I learn how to solve implicit problems in two dimensions*/
void Crank_Nicolson_scheme_Heat_2(input_2 input,int dirichlet_condition_1,int dirichlet_condition_2,int dirichlet_condition_3,double alpha)
{
    return;
}
//Backward Euler Method in one dimension (Implicit)
void BEM_Heat_1(input_1 input,int dirichlet_condition_1,int dirichlet_condition_2,int dirichlet_condition_3,double alpha){};
//Backward Euler Method in two dimensions (Implicit)
void BEM_Heat_2(input_2 input,int dirichlet_condition_1,int dirichlet_condition_2,int dirichlet_condition_3,double){};

//This will actually solve the problem and that's it.
void Heat_Eqaution_1(input_1 input)
{
    //This time maximum may be the global variable, but not all problem require us to know the time, and thus 
    //i am asking to fill the time in this function
    printf("Till what time: ");
    scanf("%lf",&time_max);
    //Now we will ask for time interval;
    printf("No of Time Steps: ");
    scanf("%d",&time_steps);

    //This will set the alpha
    double alpha;
    printf("Alpha: "); scanf("%lf",&alpha);

    /*Now is the time to take in the boudnary condition values*/
    double dirichlet_condition_1;
    double dirichlet_condition_2;
    double dirichlet_condition_3;

    printf("What do you want the Dirichlet condition assigned to number {1}: "); scanf("%lf",&dirichlet_condition_1);
    printf("What do you want the Dirichlet condition assigned to nunber {2}: "); scanf("%lf",&dirichlet_condition_2);
    printf("What do you want the Dirichlet condition assigned to nunber {3}: "); scanf("%lf",&dirichlet_condition_3);

    if(scheme_type==1)
    {
        if(scheme==1) FTCS_Heat_1(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,alpha);
        else if(scheme==2) Dufrot_Frankel_Heat_1(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,alpha);
    }
    else if(scheme_type==2)
    {
        if(scheme==1) BEM_Heat_1(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,alpha);
        else if(scheme==2) Crank_Nicolson_scheme_Heat_1(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,alpha);
    }
    return;
};

void Heat_Eqaution_2(input_2 input)
{
    printf("Till what time: "); scanf("%lf",&time_max);
    printf("No of Time Steps: "); scanf("%d",&time_steps);
    double alpha; printf("Alpha: "); scanf("%lf",&alpha);

    double dirichlet_condition_1;
    double dirichlet_condition_2;
    double dirichlet_condition_3;

    printf("What do you want the Dirichlet condition assigned to number {1}: "); scanf("%lf",&dirichlet_condition_1);
    printf("What do you want the Dirichlet condition assigned to nunber {2}: "); scanf("%lf",&dirichlet_condition_2);
    printf("What do you want the Dirichlet condition assigned to nunber {3}: "); scanf("%lf",&dirichlet_condition_3);
    if(scheme_type==1)
    {
        if(scheme==1) FTCS_Heat_2(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,alpha);
        else if(scheme==2) Dufrot_Frankel_Heat_2(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,alpha);
    }
    else if(scheme_type==2)
    {
        if(scheme==1) BEM_Heat_2(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,alpha);
        else if(scheme==2) Crank_Nicolson_scheme_Heat_2(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,alpha);
    }
    return;
}

/****************************************** LAPLACE EQUATION *****************************************/
//Explicit methods for solving laplace equations
void Explicit_sor_method_laplace_1(input_1 input,double dirchlet_condition_1,double dirchlet_condition_2,double dirchlet_condition_3)
{
    return;
}
void Explicit_sor_method_laplace_2(input_2 input,double dirchlet_condition_1,double dirchlet_condition_2,double dirchlet_condition_3)
{
    double w;
    if(debugger) printf("The scheme choosen in explicit is Successive Over Relaxation\n""If relaxation factor (w) is equal to 1 ghen then the method becomes Gauss Seidel Method\n""If 1 < w < 2 : Over Relaxation\n""If 0 < w < 1 : Under Relaxation \n");
    printf("Enter relaxation factor :");scanf("%lf",&w);
    double iterations;
    printf("What number of iteration: ");
    scanf("%lf",&iterations);

    double delta_x = L_x/grid_points_x;
    double delta_y = L_y/grid_points_y;
    double beta = delta_y/delta_x;

    //Start the clock
    clock_t start = clock();
    //Can also get the time taken for each iterations
    double* iteration_time = (double*)calloc(iterations,sizeof(double));
    clock_t i_start;
    clock_t i_end;
    
    for(int iteration = 0;iteration<iterations;iteration++)
    {
        if(computation_time) i_start = clock();
        //Now the code will traverse from the bottom left corner to the top right corner. And here it is assumed that the boundary points are given to us. So all the points at the boundary is given
        for(int j=1;j<grid_points_y-1;j++)
        {
            for(int i=1;i<grid_points_x-1;i++)
            {
                //If the interior point doesnt have any dirichlet condition then you may update it's value
                if(!input.boundary_matrix[i][j])
                {
                    input.value_matrix[i][j] = input.value_matrix[i][j] + w*((input.value_matrix[i-1][j]+input.value_matrix[i+1][j]+pow(beta,2)*(input.value_matrix[i][j+1]+input.value_matrix[i][j-1]))/(2*(1+pow(beta,2)))-input.value_matrix[i][j]);
                }
                switch (input.boundary_matrix[i][j])
                {
                case 1:
                    input.value_matrix[i][j] = dirchlet_condition_1;
                    break;
                case 2:
                    input.value_matrix[i][j] = dirchlet_condition_2;
                    break;
                case 3:
                    input.value_matrix[i][j] = dirchlet_condition_3;
                    break;
                default:
                    break;
                }
            } 
        }
        
        if(computation_time) 
        {
            i_end = clock();
            iteration_time[iteration] = (double)(i_end-i_start)/CLOCKS_PER_SEC;
        }

        //This debugger will be used to print the matrix at every iteration.
        if(debugger)
        {
            printf("Completion of iteration %d\n",iteration);
            for(int j=grid_points_y-1;j>=0;j--){
                for(int i=0;i<grid_points_x;i++)
                {
                    printf("%lf ",input.value_matrix[i][j]);
                }; 
                printf("\n");
            };
            printf("\n");
        }    
    }
    //End the Clock
    clock_t end = clock();
    if(computation_time)
    {
        double cpu_time = ((double)(end-start)/CLOCKS_PER_SEC);
        //For iteration time
        for(int iteration = 0;iteration<iterations;iteration++)
        {
            printf("Time taken for iteration %d : %f secs\n",iteration,iteration_time[iteration]);
        }
        printf("Total Time Taken : \n",cpu_time);
        free(iteration_time);
    }
    
    return;
};

// Here is the function that will be used to solve the system of linear equation using Thomas Algoritm.
//So obviously this function will be called after each row, and it will solve the equation for each row and then for each column
//This will take in two matrix. The first one being matrix A which will be 
void Implicit_sor_by_line_method_laplace_1(input_1 input,double dirchlet_condition_1,double dirchlet_condition_2,double dirchlet_condition_3)
{
    double relaxation_factor;
    printf("Relaxation factor: "); scanf("%lf",&relaxation_factor);

    return;
};

//Implicit in i
input_2 Thomas_Algorithm_sor_by_line_Method_i(input_2 input,int k,double w,double beta,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    if(debugger){
        printf("\033[1;31m");
        printf("\n");
        printf("Into the Thomas_Algorithm_Block_Iteration_Method_i\n");
        printf("input[row][%d]\n",k);
        printf("\n");
        printf("\033[0m");
    }
    //Creating Marix A
    int n = grid_points_x-2;
    /*This is a bit tricky business to deal with implicit functions and it is necessary for us to stay careful.So for that reason I ma intalising the elements with zero when forming it,. Therefore I am using calloc function rather than malloc*/
    double** A = (double**)malloc(sizeof(double*)*(n));
    //its takes in two argument, the number of elements and the size of each element. Interesting!
    for(int i=0;i<n;i++){A[i] = (double*)calloc(n,sizeof(double));}
    if(debugger){
        printf("\033[1;31mJust for debugging purpose printing the matrix A\n");printf("\033[31m");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
    }
    
    /*So here is the convention. The row is first index of A, and column is second index of A. And the top right elements will be the A[0][0], and bottom left element will be A[n-2][n-2]*/
    for(int row=0;row<n;row++){
        A[row][row] = 1;
        if(row>0) A[row][row-1] = -(w/(2*(1+pow(beta,2))));
        else if(row<n-1) A[row][row+1] = -(w/(2*(1+pow(beta,2))));
    }
    if(debugger){
        printf("\033[1;31m(Debugging)\n""Printing the matrix A\n");printf("\033[31m");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");printf("\033[0m");
        }
    }

    //Creating b matrix
    /*This is possible to create at a specific time. Since the input that we have can be used to calculate b*/
    double* b = (double*)malloc(sizeof(double)*n);
    for(int row=0;row<n;row++){b[row] = ((pow(beta,2)*w)/(2*(1+pow(beta,2))))*(input.value_matrix[k][row+2]+input.value_matrix[k][row])+(1-w)*(input.value_matrix[k][row+1]);}
    if(debugger){
        printf("\033[1;31m(Debugging)\n""Printing the matrix b after initialization\n");printf("\033[31m");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");
        printf("\033[0m");
    }
    /*Lets say that I impose the value of u_row to be eqaul to u_imposed, then I must replace A[row][row] = 1 and b[row] = u_imposed, and
    then replace all the other values of matrix A with zero in that row.*/
    for(int l=1;l<grid_points_x-1;l++)
    {
        if(input.boundary_matrix[k][l])
        {
            for(int k=0;k<n;k++) {A[l-1][k] = 0;}
            A[l-1][l-1] = 1;
            switch (input.boundary_matrix[k][l])
            {
                case 1:
                    b[l-1] = dirichlet_condition_1;
                    break;
                case 2:
                    b[l-1] = dirichlet_condition_2;
                    break;
                case 3:
                    b[l-1] = dirichlet_condition_3;
                    break;
                default:
                    break;
            }
        }
    }
    
    //Forward Elimination
    double temp;
    for(int row=0;row<n-1;row++)
    {
        //First we make the diagonal element equal to 1
        temp= A[row][row];
        for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
        b[row] = b[row]/temp;
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after making diagonal element equal to 1 at row = %d\n",row);
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after making diagonal element equal to 1 at row = %d\n",row);
            for(int i=0;i<n;i++){printf("%lf ",b[i]);}
            printf("\n");
        }
        // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
        temp = A[row+1][row];
        for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
        b[row+1] = b[row+1] - temp*b[row];
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after eliminating the element\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after eliminating the element\n");
            for(int i=0;i<n;i++){printf("%lf ",b[i]);}
            printf("\n");
        }
    }
    if(debugger){
        printf("\033[1;31m(Debugging)\n""Printing the matrix A after forward elimination\n");printf("\033[31m");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after forward elimination\n");printf("\033[31m");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");
        printf("\033[0m");
    }
    
    //Backward substitution
    input.value_matrix[k][n] = b[n-1]/A[n-1][n-1];
    for(int row=n-1;row>=1;row--)
    {
        input.value_matrix[k][row] = b[row-1]-A[row-1][row]*input.value_matrix[k][row+1];
    }
    if(debugger){
        printf("\033[1;31m(Debugging)\n""Printing the matrix A after backward substitution\n");printf("\033[31m");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[1;31m(Debugging)\n""Printing the matrix b after backward substitution\n");printf("\033[31m");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");
        printf("\033[1;31m(Debugging)\n""Printing the value matrix after backward substitution\n");printf("\033[31m");
        for(int i=0;i<grid_points_x;i++){printf("%lf ",input.value_matrix[k][i]);}
        printf("\n");
        printf("\033[0m");
    }
    //And this back substitution will finally complete all the problem, and we will be able to update all the inut values
    return input;
}
//Implicit in j
input_2 Thomas_Algorithm_sor_by_line_Method_j(input_2 input,int k,double w,double beta,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    if(debugger)
    {
        printf("\033[1;31m");
        printf("\n");
        printf("Into the Thomas_Algorithm_Block_Iteration_Method_j\n");
        printf("input[%d][column]\n",k);
        printf("\n");
        printf("\033[0m");
    }

    //Creating Marix A
    int n = grid_points_y-2;
    /*This is a bit tricky business to deal with implicit functions and it is necessary for us to stay careful.So for that reason I ma intalising the elements with zero when forming it,. Therefore I am using calloc function rather than malloc*/
    double** A = (double**)malloc(sizeof(double*)*(n));
    //its takes in two argument, the number of elements and the size of each element. Interesting!
    for(int i=0;i<n;i++){A[i] = (double*)calloc(n,sizeof(double));}
    
    if(debugger){
        printf("\033[1;31mJust for debugging purpose printing the matrix A\n");printf("\033[31m");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
    }
    
    /*So here is the convention. The row is first index of A, and column is second index of A. And the top right elements will be the A[0][0], and bottom left element will be A[n-2][n-2]*/
    for(int row=0;row<n;row++)
    {
        A[row][row] = 1;
        if(row>0) A[row][row-1] = -((w*pow(beta,2))/(2*(1+pow(beta,2))));
        else if(row<n-1) A[row][row+1] = -((w*pow(beta,2))/(2*(1+pow(beta,2))));
    }

    if(debugger){
        printf("\033[1;31m(Debugging)\n""Printing the matrix A\n");printf("\033[31m");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
    }

    //Creating b matrix
    /*This is possible to create at a specific time. Since the input that we have can be used to calculate b*/
    double* b = (double*)malloc(sizeof(double)*n);
    for(int row=0;row<n;row++){b[row] = (w/(2*(1+pow(beta,2))))*(input.value_matrix[row+2][k]+input.value_matrix[row][k])+(1-w)*(input.value_matrix[row+1][k]);}
    if(debugger){
        printf("\033[1;31m(Debugging)""Printing the matrix b after initialization\n");printf("\033[31m");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");printf("\033[0m");
    }
    /*Lets say that I impose the value of u_row to be eqaul to u_imposed, then I must replace A[row][row] = 1 and b[row] = u_imposed, and
    then replace all the other values of matrix A with zero in that row.*/
    for(int l=1;l<grid_points_y-1;l++)
    {
        if(input.boundary_matrix[l][k])
        {
            for(int c=0;c<n;c++) {A[l-1][c] = 0;}
            A[l-1][l-1] = 1;
            switch (input.boundary_matrix[l][k])
            {
                case 1:
                    b[l-1] = dirichlet_condition_1;
                    break;
                case 2:
                    b[l-1] = dirichlet_condition_2;
                    break;
                case 3:
                    b[l-1] = dirichlet_condition_3;
                    break;
                default:
                    break;
            }
        }
    }
    
    //Forward Elimination
    double temp;
    for(int row=0;row<n-1;row++)
    {
        //First we make the diagonal element equal to 1
        temp= A[row][row];
        for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
        b[row] = b[row]/temp;
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after making diagonal element equal to 1 at row = %d\n",row);
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after making diagonal element equal to 1 at row = %d\n",row);
            for(int i=0;i<n;i++){printf("%lf ",b[i]);}
            printf("\n");
        }
        // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
        temp = A[row+1][row];
        for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
        b[row+1] = b[row+1] - temp*b[row];
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after eliminating the element\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after eliminating the element\n");
            for(int i=0;i<n;i++){printf("%lf ",b[i]);}
            printf("\n");
        }
    }
    if(debugger){
        printf("\033[1;31m(Debugging)""Printing the matrix A after forward elimination\n");printf("\033[31m");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[1;31m(Debugging)\n""Printing the matrix b after forward elimination\n");printf("\033[31m");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");
        printf("\033[0m");
    }
    
    //Backward substitution
    input.value_matrix[n][k] = b[n-1]/A[n-1][n-1];
    for(int row=n-1;row>=1;row--)
    {
        input.value_matrix[row][k] = b[row-1]-A[row-1][row]*input.value_matrix[row+1][k];
    }
    if(debugger){
        printf("\033[1;31m(Debugging)\n""Printing the matrix A after backward substitution\n");printf("\033[31m\n");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                printf("%lf ",A[i][j]);
            }
            printf("\n");
        }
        printf("\033[1;31m(Debugging)\n""Printing the matrix b after backward substitution\n");printf("\033[31m\n");
        for(int i=0;i<n;i++){printf("%lf ",b[i]);}
        printf("\n");
        printf("\033[1;31m(Debugging)\n""Printing the value matrix after backward substitution\n");printf("\033[31m\n");
        for(int i=0;i<grid_points_y;i++){printf("%lf ",input.value_matrix[i][k]);}
        printf("\n");
        printf("\033[0m");
    }
    //And this back substitution will finally complete all the problem, and we will be able to update all the inut values
    return input;
}

void Implicit_sor_by_line_method_laplace_2(input_2 input,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    double relaxation_factor;
    printf("Relaxation factor: "); scanf("%lf",&relaxation_factor);
    double delta_x = L_x/grid_points_x;
    double delta_y = L_y/grid_points_y;
    double beta = delta_x/delta_y;

    double iterations;
    printf("\033[1;31m\tINSTRUCTION:\n""For one iteration two sweeps will be performed. One in I direction, and the other in J direction\n\033[0m");
    printf("Enter the iterations:"); scanf("%lf",&iterations);
    system("cls");

    clock_t start = clock();
    double* iteration_time = (double*)calloc(iterations,sizeof(double));

    for(int d = 0;d<iterations;d++)
    {
        clock_t start_i = clock();
        for(int i = 1;i<grid_points_x-1;i++)
        {
            input = Thomas_Algorithm_sor_by_line_Method_i(input,i,relaxation_factor,beta,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
            if(debugger){printf("After 'i' implicit calculation\n"); print_output_2(input);}  
        }
        for(int j = 1;j<grid_points_y-1;j++)
        {
            input = Thomas_Algorithm_sor_by_line_Method_j(input,j,relaxation_factor,beta,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
            if(debugger){printf("After 'j' implicit calculation\n"); print_output_2(input);}   
        }
        clock_t end_i = clock();
        iteration_time[d] = (double)(end_i - start_i)/CLOCKS_PER_SEC;
    }

    clock_t end = clock();
    if(computation_time)
    {
        printf("\033[1;34mTotal Time Taken = %f\n",(double)(end-start)/CLOCKS_PER_SEC);
        printf("\033[31m");
        for(int i=0;i<iterations;i++)
        {
            printf("Time Taken for iteration %d : %f\n",i,iteration_time[i]);
        }
        printf("\033[0m");
    }
    return;
}

void Laplace_Eqaution_1(input_1 input)
{
    double dirichlet_condition_1;
    double dirichlet_condition_2;
    double dirichlet_condition_3;

    printf("What do you want the Dirichlet condition assigned to number {1}: "); scanf("%lf",&dirichlet_condition_1);
    printf("What do you want the Dirichlet condition assigned to nunber {2}: "); scanf("%lf",&dirichlet_condition_2);
    printf("What do you want the Dirichlet condition assigned to nunber {3}: "); scanf("%lf",&dirichlet_condition_3);

    /*The thing with laplace equation is that it is simple. And could be reused for many kinds of problem.
    looking at the note, I see that all that i need to solve it is the intial input values and thats it. And what truly matters is the boundary conditions. */
    if(scheme_type==1)
    {
        if(scheme==1) {
            printf("You have selected \033[1m Explicit scheme\033[0m for \033[1m Laplace equation\033[0m with \033[1m Gauss-Seidel Method\033[0m\n");
            //Here will be the function call for Gauss-Seidel Method
            Explicit_sor_method_laplace_1(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
        }
    }
    else if(scheme_type==2)
    {
        printf("You have selected \033[1m Implicit scheme\033[0m for \033[1m Laplace equation\033[0m\n");
        if(scheme==1){
            printf("You have selected \033[1m Implicit scheme\033[0m for \033[1m Laplace equation\033[0m with \033[1m Block-iterative procedure\033[0m\n");
            //Here will be the function call for Block-iterative procedure
            Implicit_sor_by_line_method_laplace_1(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
        }
    }
    return;
};

void Laplace_Eqaution_2(input_2 input){
    double dirichlet_condition_1;
    double dirichlet_condition_2;
    double dirichlet_condition_3;

    printf("What do you want the Dirichlet condition assigned to number {1}: "); scanf("%lf",&dirichlet_condition_1);
    printf("What do you want the Dirichlet condition assigned to nunber {2}: "); scanf("%lf",&dirichlet_condition_2);
    printf("What do you want the Dirichlet condition assigned to nunber {3}: "); scanf("%lf",&dirichlet_condition_3);

    /*The thing with laplace equation is that it is simple. And could be reused for many kinds of problem.
    looking at the note, I see that all that i need to solve it is the intial input values and thats it. And what truly matters is the boundary conditions. */
    if(scheme_type==1)
    {
        if(scheme==1) {
            printf("You have selected \033[1m Explicit scheme\033[0m for \033[1m Laplace equation\033[0m with \033[1m Gauss-Seidel Method\033[0m\n");
            //Here will be the function call for Gauss-Seidel Method
            Explicit_sor_method_laplace_2(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
        }
    }
    else if(scheme_type==2)
    {
        printf("You have selected \033[1m Implicit scheme\033[0m for \033[1m Laplace equation\033[0m\n");
        if(scheme==1){
            printf("You have selected \033[1m Implicit scheme\033[0m for \033[1m Laplace equation\033[0m with \033[1m Block-iterative procedure\033[0m\n");
            //Here will be the function call for Block-iterative procedure
            Implicit_sor_by_line_method_laplace_2(input,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
        }
    }
    return;
};

/******************************** WAVE EQUATION ****************************** */

/*This section include algorithms for solving the wave equation, both linear and non linear,
And the problem will only be solved for one dimensional*/

//This includes Explicit Methods for linear

//This will be the Euler FTFS mehtod for solving wave equation in 1D explicitly.
double Euler_FTFS_explicit_wave_1(input_1 input,int i,double corant_number){
    return input.value_matrix[i] - (corant_number)*(input.value_matrix[i+1] - input.value_matrix[i]);
}
//This will be the fucntion for Euler FTCS method
double Euler_FTCS_explicit_Wave_1(input_1 input,int i,double corant_number){
    return input.value_matrix[i] - (corant_number/2)*(input.value_matrix[i+1] + input.value_matrix[i-1]);
}
//The first upwinding scheme
double Euler_first_upwinding_explicit_wave_1(input_1 input, int i, double corant_number){
    if(corant_number>0) return (input.value_matrix[i] - corant_number*(input.value_matrix[i] - input.value_matrix[i-1]));
    else if(corant_number<0) return input.value_matrix[i] - corant_number*(input.value_matrix[i+1] - input.value_matrix[i]);
    return input.value_matrix[i];
}
//The lax method
double Lax_explicit_wave_1(input_1 input,int i,double corant_number)
{
    return (0.5*(input.value_matrix[i+1] + input.value_matrix[i-1]) - (corant_number/2)*(input.value_matrix[i+1] - input.value_matrix[i-1]));
}
//The lax wendroff method
double Lax_wendroff_explicit_wave_1(input_1 input,int i, double corant_number)
{
    return input.value_matrix[i]-(corant_number/2)*(input.value_matrix[i+1]- input.value_matrix[i-1]) + (0.5*pow(corant_number,2))*(input.value_matrix[i+1]+input.value_matrix[i-1]-2*input.value_matrix[i]);
}

input_1 Euler_Explicit_Linear_Wave_1(input_1 input,double corant_number,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3,double f(input_1,int,double))
{
    if(debugger) printf("\033[1;31mInto the function Euler_FTFS_explicit_wave_1\033[0m\n");
    for(int t=0;t<time_steps;t++)
    {
        //if(file_output) fprintf(fr,"%lf,",input.value_matrix[0]);
        if(debugger) printf("Time step = %d\n",t+1);
        for(int i=1;i<grid_points_x-1;i++)
        {
            //if you see a boundary point assign the dirichlet condition
            if(input.boundary_matrix[i]!=0)
            {
                switch (input.boundary_matrix[i])
                {
                case 1:
                    input.value_matrix[i] = dirichlet_condition_1;
                    break;
                case 2:
                    input.value_matrix[i] = dirichlet_condition_2;
                case 3:
                    input.value_matrix[i] = dirichlet_condition_3;
                default:
                    break;
                }
            }
            //If not update the point according to scheme
            else{input.value_matrix[i] = f(input,i,corant_number);}
            if(debugger) printf("input.value_matrix[%d] %lf\n",i,input.value_matrix[i]);
            //if(file_output) { fprintf(fr,"%lf,",input.value_matrix[i]);}
        }
        //if(file_output) fprintf(fr,"%lf\n",input.value_matrix[grid_points_x-1]);
    }
    return input;
}

input_1 Lax_Wendroff_Explicit_Multistep_Linear_Wave_1(input_1 input, double corant_number,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    printf("\033[31mThe Corant number for stability must be less than 1\033[0m\n"
        "Given courant number :\033[1;31m %lf\033[0m\n",corant_number);
    if(corant_number>1)
    {
        printf("\033[31mThe corant number that is entered is more than 1, and this may lead to unstable results.\033[0m\n");
    }
    double* intermediate_value_matrix = (double*)calloc(grid_points_x,sizeof(double));
    for(int t=0;t<time_steps;t++)
    {
        //if(file_output) fprintf(fr,"%lf,",input.value_matrix[0]);
        if(debugger) printf("Time step = %d\n",t+1);
        //This will update the intermediate values of the matrix
        for(int i=1;i<grid_points_x-1;i++){
            intermediate_value_matrix[i] = 0.5*(input.value_matrix[i+1]+input.value_matrix[i]) - (corant_number/2)*(input.value_matrix[i+1]- input.value_matrix[i]);
        }
        //This will update the next time step
        for(int i=1;i<grid_points_x-1;i++){
            if(input.boundary_matrix[i]!=0)
            {
                switch (input.boundary_matrix[i])
                {
                case 1:
                    input.value_matrix[i] = dirichlet_condition_1;
                    break;
                case 2:
                    input.value_matrix[i] = dirichlet_condition_2;
                    break;
                case 3:
                    input.value_matrix[i] = dirichlet_condition_3;
                    break;
                default:
                    break;
                }
            }
            else 
            {
                input.value_matrix[i] = input.value_matrix[i] - corant_number*(intermediate_value_matrix[i+1]-intermediate_value_matrix[i-1]);
            }
            if(debugger) printf("input.value_matrix[%d] %lf\n",i,input.value_matrix[i]);
            //if(file_output) { fprintf(fr,"%lf,",input.value_matrix[i]);}
        }
        //if(file_output) fprintf(fr,"%lf\n",input.value_matrix[grid_points_x-1]);
    }
    free(intermediate_value_matrix);
    return input;
}

input_1 Mac_Cormack_Explicit_method_Linear_Wave_1(input_1 input, double corant_number,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    printf("\033[31mThe Corant number for stability must be less than 1\033[0m\n"
        "Given courant number :\033[1;31m %lf\033[0m\n",corant_number);
    if(corant_number>1)
    {
        printf("\033[31mThe corant number that is entered is more than 1, and this may lead to unstable results.\033[0m\n");
    }
    double* intermediate_value_matrix = (double*)calloc(grid_points_x,sizeof(double));
    for(int t=0;t<time_steps;t++)
    {
        //if(file_output) fprintf(fr,"%lf,",input.value_matrix[0]);
        if(debugger) printf("Time step = %d\n",t+1);
        //This will update the intermediate values of the matrix
        for(int i=1;i<grid_points_x-1;i++){
            intermediate_value_matrix[i] = input.value_matrix[i] - (corant_number)*(input.value_matrix[i+1]-input.value_matrix[i]);
        }
        //This will update the next time step
        for(int i=1;i<grid_points_x-1;i++){
            if(input.boundary_matrix[i]!=0)
            {
                switch (input.boundary_matrix[i])
                {
                case 1:
                    input.value_matrix[i] = dirichlet_condition_1;
                    break;
                case 2:
                    input.value_matrix[i] = dirichlet_condition_2;
                    break;
                case 3:
                    input.value_matrix[i] = dirichlet_condition_3;
                    break;
                default:
                    break;
                }
            }
            else 
            {
                input.value_matrix[i] = 0.5*((input.value_matrix[i]+intermediate_value_matrix[i]) - corant_number*(intermediate_value_matrix[i]-intermediate_value_matrix[i-1]));
            }
            if(debugger) printf("input.value_matrix[%d] %lf\n",i,input.value_matrix[i]);
            //if(file_output) { fprintf(fr,"%lf,",input.value_matrix[i]);}
        }
        //if(file_output) fprintf(fr,"%lf\n",input.value_matrix[grid_points_x-1]);
    }
    free(intermediate_value_matrix);
    return input;
}

//This includes Implicit Methods for linear
double** A_matrix_Crank_Nicolsen_Implicit(double **A,double corant_number)
{
    for(int row=0;row<grid_points_x-2;row++)
    {
        A[row][row] = 1;
        if(row==0)  A[row][row+1] = (0.25)*(corant_number);
        else if(row==grid_points_x-3) A[row][row-1] = (-0.25)*(corant_number);
        else{
            A[row][row+1] = (0.25)*(corant_number);
            A[row][row-1] = (-0.25)*(corant_number);
        }
    }
    return A;
}
double* B_matrix_Crank_Nicolsen_Implicit(input_1 input,double *B,double corant_number)
{
    for(int row=0;row<grid_points_x-2;row++)
    {
        if(row==0)  B[row] = input.value_matrix[row+1]-(corant_number/4)*(input.value_matrix[row+2]-input.value_matrix[row])+(corant_number/4)*input.value_matrix[row];
        else if(row==grid_points_x-3) B[row] = input.value_matrix[row+1]-(corant_number/4)*(input.value_matrix[row+2]-input.value_matrix[row])-(corant_number/4)*input.value_matrix[row+2];
        else B[row] = input.value_matrix[row+1]-(corant_number/4)*(input.value_matrix[row+2]-input.value_matrix[row]);
    }
    return B;
}
double** A_matrix_Euler_BTCS_Implicit(double **A,double corant_number)
{
    for(int row=0;row<grid_points_x-2;row++)
    {
        A[row][row] = 1;
        if(row==0)  A[row][row+1] = (0.5)*(corant_number);
        else if(row==grid_points_x-3) A[row][row-1] = (-0.5)*(corant_number);
        else{
             A[row][row+1] = (0.5)*(corant_number);
             A[row][row-1] = (-0.5)*(corant_number);
        }
    }
    return A;
}
double* B_matrix_Euler_BTCS_Implicit(input_1 input,double *B,double corant_number)
{
    for(int row=0;row<grid_points_x-2;row++)
    {
        if(row == 0) B[row] = input.value_matrix[row+1] - (-corant_number/2)*input.value_matrix[row];
        else if(row == grid_points_x-3) B[row] = input.value_matrix[row+1] - (corant_number/2)*input.value_matrix[row+2];
        else{
            B[row] = input.value_matrix[row+1];
        }
    }
    return B;
}

void Euler_BTCS_Implicit_Linear_Wave_1(input_1 input,double corant_number,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    const int n = grid_points_x - 2;
    for(int t=0;t<time_steps;t++)
    {
        if(debugger){printf("Time step = %d\n",t);}
        //construction of A matrix
        double** A = (double**)calloc(n,sizeof(double*));
        for(int i=0;i<n;i++) A[i] = (double*)calloc(n,sizeof(double));
        A = A_matrix_Euler_BTCS_Implicit(A,corant_number);

        if(debugger){
            printf("Just for debugging purpose printing the matrix A\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
        }
        
        // Construction of b
        double* B = (double*)calloc(n,sizeof(double));
        B = B_matrix_Euler_BTCS_Implicit(input,B,corant_number);
        
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after initialization\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
        }
        
        //Imposing boundary condtion
        for(int l=1;l<grid_points_x-1;l++)
        {
            if(input.boundary_matrix[l])
            {
                for(int k=0;k<n;k++) {A[l-1][k] = 0;}
                A[l-1][l-1] = 1;
                switch (input.boundary_matrix[l])
                {
                    case 1:
                        B[l-1] = dirichlet_condition_1;
                        break;
                    case 2:
                        B[l-1] = dirichlet_condition_2;
                        break;
                    case 3:
                        B[l-1] = dirichlet_condition_3;
                        break;
                    default:
                        break;
                }
            }
        } 
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after imposing boundary condition\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after imposing boundary condition\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
        }
        
        //Forward Elimination
        double temp;
        for(int row=0;row<n-1;row++)
        {
            //First we make the diagonal element equal to 1
            temp= A[row][row];
            for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
            B[row] = B[row]/temp;
            if(debugger){
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after making diagonal element equal to 1 at row = %d\n",row);
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        printf("%lf ",A[i][j]);
                    }
                    printf("\n");
                }
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after making diagonal element equal to 1 at row = %d\n",row);
                for(int i=0;i<n;i++){printf("%lf ",B[i]);}
                printf("\n");
            }
            // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
            temp = A[row+1][row];
            for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
            B[row+1] = B[row+1] - temp*B[row];
            if(debugger){
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after eliminating the element\n");
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        printf("%lf ",A[i][j]);
                    }
                    printf("\n");
                }
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after eliminating the element\n");
                for(int i=0;i<n;i++){printf("%lf ",B[i]);}
                printf("\n");
            }
        }
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after forward elimination\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after forward elimination\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
        }
        
        //Backward substitution
        input.value_matrix[n] = B[n-1]/A[n-1][n-1];
        for(int row=n-1;row>=1;row--)
        {
            input.value_matrix[row] = B[row-1]-A[row-1][row]*input.value_matrix[row+1];
        }
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after backward substitution\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after backward substitution\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the value matrix after backward substitution\n");
            for(int i=0;i<grid_points_x;i++){printf("%lf ",input.value_matrix[i]);}
            printf("\n");
        }
        
        if(file_output)
        {
            for(int i=0;i<grid_points_x;i++){
                //fprintf(fr,"%lf,",input.value_matrix[i]);
                //if(i == grid_points_x-1) fprintf(fr,"%lf\n",input.value_matrix[i]);
            }
        }
        free(A);
        free(B);
    }
    return;
}
void Crank_Nicolsen_Implicit_Linear_Waver_1(input_1 input,double corant_number,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    int n = grid_points_x - 2;
    for(int t=0;t<time_steps;t++)
    {
        if(debugger){printf("Time step = %d\n",t);}
        //construction of A matrix
        double** A = (double**)calloc(n,sizeof(double*));
        for(int i=0;i<n;i++) A[i] = (double*)calloc(n,sizeof(double));
        A = A_matrix_Crank_Nicolsen_Implicit(A,corant_number);
        if(debugger){
            printf("Just for debugging purpose printing the matrix A\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
        }
        
        // Construction of b
        double* B = (double*)calloc(n,sizeof(double));
        B = B_matrix_Crank_Nicolsen_Implicit(input,B,corant_number);
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after initialization\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
        }
        
        //Imposing boundary condtion
        for(int l=1;l<grid_points_x-1;l++)
        {
            if(input.boundary_matrix[l])
            {
                for(int k=0;k<n;k++) {A[l-1][k] = 0;}
                A[l-1][l-1] = 1;
                switch (input.boundary_matrix[l])
                {
                    case 1:
                        B[l-1] = dirichlet_condition_1;
                        break;
                    case 2:
                        B[l-1] = dirichlet_condition_2;
                        break;
                    case 3:
                        B[l-1] = dirichlet_condition_3;
                        break;
                    default:
                        break;
                }
            }
        } 
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after imposing boundary condition\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after imposing boundary condition\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
        }
        
        //Forward Elimination
        double temp;
        for(int row=0;row<n-1;row++)
        {
            //First we make the diagonal element equal to 1
            temp= A[row][row];
            for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
            B[row] = B[row]/temp;
            if(debugger){
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after making diagonal element equal to 1 at row = %d\n",row);
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        printf("%lf ",A[i][j]);
                    }
                    printf("\n");
                }
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after making diagonal element equal to 1 at row = %d\n",row);
                for(int i=0;i<n;i++){printf("%lf ",B[i]);}
                printf("\n");
            }
            // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
            temp = A[row+1][row];
            for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
            B[row+1] = B[row+1] - temp*B[row];
            if(debugger){
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after eliminating the element\n");
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        printf("%lf ",A[i][j]);
                    }
                    printf("\n");
                }
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after eliminating the element\n");
                for(int i=0;i<n;i++){printf("%lf ",B[i]);}
                printf("\n");
            }
        }
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after forward elimination\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after forward elimination\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
        }
        
        //Backward substitution
        input.value_matrix[n] = B[n-1]/A[n-1][n-1];
        for(int row=n-1;row>=1;row--)
        {
            input.value_matrix[row] = B[row-1]-A[row-1][row]*input.value_matrix[row+1];
        }
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after backward substitution\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after backward substitution\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the value matrix after backward substitution\n");
            for(int i=0;i<grid_points_x;i++){printf("%lf ",input.value_matrix[i]);}
            printf("\n");
        }
        
        if(file_output)
        {
            for(int i=0;i<grid_points_x;i++){
                //fprintf(fr,"%lf,",input.value_matrix[i]);
                //if(i == grid_points_x-1) fprintf(fr,"%lf\n",input.value_matrix[i]);
            }
        }
        
        free(A);
        free(B);
    }
    return;
}


//This includes Explicit Methods for Non Linear
void Lax_Explicit_Non_Linear_Wave_1(input_1 input,double corant_number,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    if(debugger) printf("\033[1;31mInto the function Lax_Explicit_Non_Linear_Wave_1\033[0m\n");
    for(int t=0;t<time_steps;t++)
    {
        //if(file_output) fprintf(fr,"%lf,",input.value_matrix[0]);
        if(debugger) printf("Time step = %d\n",t+1);
        for(int i=1;i<grid_points_x-1;i++)
        {
            //if you see a boundary point assign the dirichlet condition
            if(input.boundary_matrix[i]!=0)
            {
                switch (input.boundary_matrix[i])
                {
                case 1:
                    input.value_matrix[i] = dirichlet_condition_1;
                    break;
                case 2:
                    input.value_matrix[i] = dirichlet_condition_2;
                    break;
                case 3:
                    input.value_matrix[i] = dirichlet_condition_3;
                    break;
                default:
                    break;
                }
            }
            //If not update the point according to scheme
            else{
                input.value_matrix[i] = 0.5*(input.value_matrix[i+1]+input.value_matrix[i-1])-(corant_number/4)*(pow(input.value_matrix[i+1],2) - pow(input.value_matrix[i-1],2));
            }
            if(debugger) printf("input.value_matrix[%d] %lf\n",i,input.value_matrix[i]);
            //if(file_output) { fprintf(fr,"%lf,",input.value_matrix[i]);}
        }
        //if(file_output) fprintf(fr,"%lf\n",input.value_matrix[grid_points_x-1]);
    }
    return;
}
void Lax_Wendroff_Explicit_Non_Linear_Wave_1(input_1 input,double corant_number,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    if(debugger) printf("\033[1;31mInto the function Euler_FTFS_explicit_wave_1\033[0m\n");
    for(int t=0;t<time_steps;t++)
    {
        //if(file_output) fprintf(fr,"%lf,",input.value_matrix[0]);
        if(debugger) printf("Time step = %d\n",t+1);
        for(int i=1;i<grid_points_x-1;i++)
        {
            //if you see a boundary point assign the dirichlet condition
            if(input.boundary_matrix[i]!=0)
            {
                switch (input.boundary_matrix[i])
                {
                case 1:
                    input.value_matrix[i] = dirichlet_condition_1;
                    break;
                case 2:
                    input.value_matrix[i] = dirichlet_condition_2;
                case 3:
                    input.value_matrix[i] = dirichlet_condition_3;
                default:
                    break;
                }
            }
            //If not update the point according to scheme
            else{
                double first_term = (-corant_number/2)*(pow(input.value_matrix[i+1],2) - pow(input.value_matrix[i-1],2));
                double second_term = (pow(corant_number,2)/4)*((input.value_matrix[i+1]-input.value_matrix[i])*(pow(input.value_matrix[i+1],2) - pow(input.value_matrix[i],2)));
                double third_term = (-pow(corant_number,2)/4)*((input.value_matrix[i]-input.value_matrix[i-1])*(pow(input.value_matrix[i],2) - pow(input.value_matrix[i-1],2)));
                input.value_matrix[i] = input.value_matrix[i] + first_term + second_term + third_term;
            }
            if(debugger) printf("input.value_matrix[%d] %lf\n",i,input.value_matrix[i]);
            //if(file_output) { fprintf(fr,"%lf,",input.value_matrix[i]);}
        }
        //if(file_output) fprintf(fr,"%lf\n",input.value_matrix[grid_points_x-1]);
    }
    return;
}
void Mac_Cormack_Explicit_Non_Linear_Wave_1(input_1 input,double corant_number,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    if(debugger) printf("\033[1;31mInto the function Mac_Cormack_Explicit_Non_Linear_Wave_1\033[0m\n");
    for(int t=0;t<time_steps;t++)
    {
        //if(file_output) fprintf(fr,"%lf,",input.value_matrix[0]);
        if(debugger) printf("Time step = %d\n",t+1);
        double* intermediate_value_matrix = (double*)calloc(grid_points_x,sizeof(double));
        if(debugger) printf("Time Step = %d\n",t);
        //This will update the intermediate values of the matrix
        for(int i=1;i<grid_points_x-1;i++){
            intermediate_value_matrix[i] = input.value_matrix[i] - (corant_number)*((pow(input.value_matrix[i+1],2)/2)-(pow(input.value_matrix[i],2)/2));
        }
        //This will update the next time step
        for(int i=1;i<grid_points_x-1;i++){
            if(input.boundary_matrix[i]!=0)
            {
                switch (input.boundary_matrix[i])
                {
                case 1:
                    input.value_matrix[i] = dirichlet_condition_1;
                    break;
                case 2:
                    input.value_matrix[i] = dirichlet_condition_2;
                    break;
                case 3:
                    input.value_matrix[i] = dirichlet_condition_3;
                    break;
                default:
                    break;
                }
            }
            else 
            {input.value_matrix[i] = 0.5*((input.value_matrix[i]+intermediate_value_matrix[i]) - corant_number*((pow(input.value_matrix[i],2)/2)-(pow(input.value_matrix[i-1],2)/2)));}
            if(debugger) printf("input.value_matrix[%d] %lf\n",i,input.value_matrix[i]);
            //if(file_output) { fprintf(fr,"%lf,",input.value_matrix[i]);}
        }
        //if(file_output) fprintf(fr,"%lf\n",input.value_matrix[grid_points_x-1]);
        free(intermediate_value_matrix);
    }
    return;
}


double** A_matrix_Beam_and_Warming_Implicit(input_1 input,double **A,double corant_number){
    for(int row=0;row<grid_points_x-2;row++)
    {
        A[row][row] = 1;
        if(row==0)  A[row][row+1] = (-0.25*corant_number)*input.value_matrix[row];
        else if(row==grid_points_x-3) A[row][row-1] = (0.25*corant_number)*input.value_matrix[row+2];
        else{
            A[row][row+1] = (-0.25*corant_number)*input.value_matrix[row];
            A[row][row-1] = (0.25*corant_number)*input.value_matrix[row+2];
        }
    }
    return A;
}
double* B_matrix_Beam_and_Warming_Implicit(input_1 input,double *B,double corant_number){
    for(int row=0;row<grid_points_x-2;row++)
    {
        B[row] = input.value_matrix[row+1] - (corant_number/2)*(pow(input.value_matrix[row+2],2) - pow(input.value_matrix[row],2)) + (corant_number/4)*pow(input.value_matrix[row+2],2) - (corant_number/4)*pow(input.value_matrix[row],2);
        if(row == 0) B[row] = B[row] - (-0.25*corant_number)*input.value_matrix[row]*input.value_matrix[row];
        else if(row == grid_points_x-3) B[row] = B[row] - (0.25*corant_number)*input.value_matrix[row+2]*input.value_matrix[row];
    }
    return B;
}
//This includes Implicit methods for Non Linear
void Beam_and_Warming_Implicit_Non_Linear_Wave_1(input_1 input,double corant_number,double dirichlet_condition_1,double dirichlet_condition_2,double dirichlet_condition_3)
{
    const int n = grid_points_x - 2;
    for(int t=0;t<time_steps;t++)
    {
        if(debugger){printf("Time step = %d\n",t);}
        //construction of A matrix
        double** A = (double**)calloc(n,sizeof(double*));
        for(int i=0;i<n;i++) A[i] = (double*)calloc(n,sizeof(double));
        A = A_matrix_Beam_and_Warming_Implicit(input,A,corant_number);

        if(debugger){
            printf("Just for debugging purpose printing the matrix A\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
        }
        
        // Construction of b
        double* B = (double*)calloc(n,sizeof(double));
        B = B_matrix_Beam_and_Warming_Implicit(input,B,corant_number);
        
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after initialization\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
        }
        
        //Imposing boundary condtion
        for(int l=1;l<grid_points_x-1;l++)
        {
            if(input.boundary_matrix[l])
            {
                for(int k=0;k<n;k++) {A[l-1][k] = 0;}
                A[l-1][l-1] = 1;
                switch (input.boundary_matrix[l])
                {
                    case 1:
                        B[l-1] = dirichlet_condition_1;
                        break;
                    case 2:
                        B[l-1] = dirichlet_condition_2;
                        break;
                    case 3:
                        B[l-1] = dirichlet_condition_3;
                        break;
                    default:
                        break;
                }
            }
        } 
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after imposing boundary condition\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after imposing boundary condition\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
        }
        
        //Forward Elimination
        double temp;
        for(int row=0;row<n-1;row++)
        {
            //First we make the diagonal element equal to 1
            temp= A[row][row];
            for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
            B[row] = B[row]/temp;
            if(debugger){
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after making diagonal element equal to 1 at row = %d\n",row);
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        printf("%lf ",A[i][j]);
                    }
                    printf("\n");
                }
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after making diagonal element equal to 1 at row = %d\n",row);
                for(int i=0;i<n;i++){printf("%lf ",B[i]);}
                printf("\n");
            }
            // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
            temp = A[row+1][row];
            for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
            B[row+1] = B[row+1] - temp*B[row];
            if(debugger){
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after eliminating the element\n");
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        printf("%lf ",A[i][j]);
                    }
                    printf("\n");
                }
                printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after eliminating the element\n");
                for(int i=0;i<n;i++){printf("%lf ",B[i]);}
                printf("\n");
            }
        }
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after forward elimination\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after forward elimination\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
        }
        
        //Backward substitution
        input.value_matrix[n] = B[n-1]/A[n-1][n-1];
        for(int row=n-1;row>=1;row--)
        {
            input.value_matrix[row] = B[row-1]-A[row-1][row]*input.value_matrix[row+1];
        }
        if(debugger){
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix A after backward substitution\n");
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    printf("%lf ",A[i][j]);
                }
                printf("\n");
            }
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the matrix b after backward substitution\n");
            for(int i=0;i<n;i++){printf("%lf ",B[i]);}
            printf("\n");
            printf("\033[1;31m(Debugging)\033[0m\n""Printing the value matrix after backward substitution\n");
            for(int i=0;i<grid_points_x;i++){printf("%lf ",input.value_matrix[i]);}
            printf("\n");
        }
        
        if(file_output)
        {
            for(int i=0;i<grid_points_x;i++){
                //fprintf(fr,"%lf,",input.value_matrix[i]);
                //if(i == grid_points_x-1) fprintf(fr,"%lf\n",input.value_matrix[i]);
            }
        }
        free(A);
        free(B);
    }
    return;
}

input_1 Wave_Eqaution_1(input_1 input, int islinear)
{
    printf("Into to the Wave Equation 1 Function\n");

    printf("Enter Total Time: "); scanf("%lf",&time_max);
    printf("Enter Time steps: "); scanf("%d",&time_steps);
    double delta_x = L_x/grid_points_x;
    double delta_t = time_max/time_steps;
    
    // Boundary condition
    double dirichlet_condition_1;
    double dirichlet_condition_2;
    double dirichlet_condition_3;
    printf("What do you want the Dirichlet condition assigned to number {1}: "); scanf("%lf",&dirichlet_condition_1);
    printf("What do you want the Dirichlet condition assigned to nunber {2}: "); scanf("%lf",&dirichlet_condition_2);
    printf("What do you want the Dirichlet condition assigned to nunber {3}: "); scanf("%lf",&dirichlet_condition_3);

    
    if(islinear)
    {
        // Corant number and value of wave speed (a) exist only in linear wave eqaution.
        double a;
        printf("Enter a: "); scanf("%lf",&a);
        printf("The value entered of a is : \n",a);
        double corant_number = (a*delta_t)/(delta_x);

        printf("You have selected Linear Wave Eqaution\n");
        //if the scheme is Explicit
        if(scheme_type == 1)
        {
            printf("You have selected \033[1;32mExplicit Scheme\033[0m\n");
            if(scheme == 1){
                printf("Selected \033[1;32mEuler FTFS method\033[1;32m\n");
                input = Euler_Explicit_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,Euler_FTFS_explicit_wave_1);
            }
            else if(scheme == 2){
                printf("Selected \033[1;32mEuler FTCS method\033[1;32m\n");
                input = Euler_Explicit_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,Euler_FTCS_explicit_Wave_1);
            }
            else if(scheme == 3){
                printf("Selected \033[1;32mThe First Upwinding Differencing method\033[1;32m\n");
                input = Euler_Explicit_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,Euler_first_upwinding_explicit_wave_1);
            }
            else if(scheme == 4){
                printf("Selected \033[1;32mLax method\033[1;32m\n");
                input = Euler_Explicit_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,Lax_explicit_wave_1);
            }
            else if(scheme == 5){
                printf("Selected \033[1;32mLax Wendroff method\033[1;32m\n");
                input = Euler_Explicit_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3,Lax_wendroff_explicit_wave_1);
            }
            else if(scheme == 6){
                printf("Selected \033[1;32mLax Wendroff MultiStep method\033[1;32m\n");
                input = Lax_Wendroff_Explicit_Multistep_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
            }
            else if(scheme == 7)
            {
                printf("Selected \033[1;32mMac Cormack method\033[1;32m\n");
                input = Mac_Cormack_Explicit_method_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
            }
        }
        //If the scheme is Implicit
        else if(scheme_type == 2)
        {
            printf("You have selected \033[1;32mImplicit Scheme\033[0m\n");
            switch (scheme)
            {
                case 1:
                    printf("Selected \033[1;32mEuler's BTCS Method\0333[1;32m\n");
                    Euler_BTCS_Implicit_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
                    break;
                case 2:
                    printf("Selected \033[1;32mCrank Nicholsen Method\0333[1;32m\n");
                    Crank_Nicolsen_Implicit_Linear_Waver_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
                    break;
                default:
                    break;
            }
        }
    } else {
        double corant_number = (delta_t)/(delta_x);
        printf("You have selected Non linear Wave Eqaution\n");
        //If Explicit scheme
        if(scheme_type == 1)
        {
            printf("You have selected \033[1;32mExplicit Scheme\033[0m\n");
            switch (scheme)
            {
            case 1:
                printf("Selected \033[1;32mLax Method\033[1;32m\n");
                Lax_Explicit_Non_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
                break;
            case 2:
                printf("Selected \033[1;32mLax Wendroff Method\033[1;32m\n");
                Lax_Wendroff_Explicit_Non_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
                break;
            case 3:
                printf("Selected \033[1;32mMac Cormack Method\033[1;32m\n");
                Mac_Cormack_Explicit_Non_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
                break;
            default:
                break;
            }
        }
        //if Implicit Scheme
        else if(scheme_type == 2)
        {
            printf("You have selected \033[1;32mImplicit Scheme\033[0m\n");
            switch (scheme)
            {
                case 1:
                    printf("Selected \033[1;32mBeam And Warning Implicit Method\033[1;32m\n");
                    Beam_and_Warming_Implicit_Non_Linear_Wave_1(input,corant_number,dirichlet_condition_1,dirichlet_condition_2,dirichlet_condition_3);
                    break;
                default:
                    break;
            }
        }
    }
    return input;
}

/* ANIMATION ********************************************************************************** */

void opening_animation_2()
{
    printf("\033[?25l");//Hide the curser
    printf("\033[H"); //Clear
    //This section is for the static picture of the OpenFOAM type interface.
    {
        //This will be used to generate the top edge of the box.
        printf("/*");
        for(int i=0;i<60;i++){printf("-");}
        printf("*\\");
        printf("\n");

        //Now we will add the equal to sign for the top if the triangle
        printf("| ");
        for(int equal=0;equal<10;equal++) printf("=");
        for(int equal=0;equal<14;equal++) printf(" ");
        printf(" |");
        printf("                                   |");
        printf("\n");

        //There are a total of 5 lines (rows) and this will be noted with the variable row
        for(int row=0;row<4;row++)
        {
            printf("|");
            //Now is the first row we will decided how many spaces must be added before printing the traingle
            for(int space=0;space<row+1;space++){printf(" ");}
            //And then we will print the double edge of the trianlge
            printf("\\\\");
            //then depending on the row, we will again add the space
            for(int space=0;space<2*(4-row)-1;space++){printf(" ");}
            //And then we print the other edge of the trinlge
            printf("/");
            //After the trianlge is printed we will again print the space
            for(int space=0;space<row+1;space++){printf(" ");}
            //Now we type the open form full form
            switch (row)
            {
            case 0:
                printf("F ield        ");
                //No we may add the closing sign
                break;
            case 1:
                printf("O peration    ");
                break;
            case 2:
                printf("A nd          ");
                break;
            case 3:
                printf("M anipulation ");
                break;
            default:
                break;
            }
            printf("|");
            printf("                                   |");
            printf("\n");
        };
        
        //And finally the bottom edge of the box
        //This will be used to generate the top edge of the box.
        printf("\\*");
        for(int i=0;i<60;i++){printf("-");}
        printf("*/");
        printf("\n");
        Sleep(2000);
    }   
    system("cls");
    //This section is for the bottom edge to disapper
    {
        for(int bottom_row=0;bottom_row<32;bottom_row++)
        {
            printf("\033[H");
            //This will be used to generate the top edge of the box.
            printf("/*");
            for(int i=0;i<60;i++){printf("-");}
            printf("*\\");
            printf("\n");

            //Now we will add the equal to sign for the top if the triangle
            printf("| ");
            for(int equal=0;equal<10;equal++) printf("=");
            for(int equal=0;equal<14;equal++) printf(" ");
            printf(" |");
            printf("                                   |");
            printf("\n");

            //There are a total of 5 lines (rows) and this will be noted with the variable row
            for(int row=0;row<4;row++)
            {
                printf("|");
                //Now is the first row we will decided how many spaces must be added before printing the traingle
                for(int space=0;space<row+1;space++){printf(" ");}
                //And then we will print the double edge of the trianlge
                printf("\\\\");
                //then depending on the row, we will again add the space
                for(int space=0;space<2*(4-row)-1;space++){printf(" ");}
                //And then we print the other edge of the trinlge
                printf("/");
                //After the trianlge is printed we will again print the space
                for(int space=0;space<row+1;space++){printf(" ");}
                //Now we type the open form full form
                switch (row)
                {
                case 0:
                    printf("F ield        ");
                    //Now we may add the closing sign
                    break;
                case 1:
                    printf("O peration    ");
                    break;
                case 2:
                    printf("A nd          ");
                    break;
                case 3:
                    printf("M anipulation ");
                    break;
                default:
                    break;
                }
                printf("|");
                printf("                                   |");
                printf("\n");
            };
            
            //And finally the bottom edge of the box
            //This will be used to generate the top edge of the box.
            printf("\\*");
            for(int i=0;i<(30-bottom_row);i++){printf("-");}
            for(int i=0;i<2*bottom_row;i++){printf(" ");}
            for(int i=0;i<(30-bottom_row);i++){printf("-");}
            printf("*/");
            printf("\n");
            Sleep(25);
        }
    }
    system("cls");
    //This section is for the side edge to disapper
    {
        for(int side_edge=0;side_edge<5;side_edge++){
            printf("\033[H");
            //This will be used to generate the top edge of the box.
            printf("/*");
            for(int i=0;i<60;i++){printf("-");}
            printf("*\\");
            printf("\n");

            //Now we will add the equal to sign for the top if the triangle
            printf("| ");
            for(int equal=0;equal<10;equal++) printf("=");
            for(int equal=0;equal<14;equal++) printf(" ");
            printf(" |");
            printf("                                   |");
            printf("\n");

            //There are a total of 5 lines (rows) and this will be noted with the variable row
            for(int row=0;row<4;row++)
            {
                if((row+side_edge)<4) printf("|");
                else printf(" ");
                //Now is the first row we will decided how many spaces must be added before printing the traingle
                for(int space=0;space<row+1;space++){printf(" ");}
                //And then we will print the double edge of the trianlge
                printf("\\\\");
                //then depending on the row, we will again add the space
                for(int space=0;space<2*(4-row)-1;space++){printf(" ");}
                //And then we print the other edge of the trinlge
                printf("/");
                //After the trianlge is printed we will again print the space
                for(int space=0;space<row+1;space++){printf(" ");}
                //Now we type the open form full form
                switch (row)
                {
                case 0:
                    printf("F ield        ");
                    //Now we may add the closing sign
                    break;
                case 1:
                    printf("O peration    ");
                    break;
                case 2:
                    printf("A nd          ");
                    break;
                case 3:
                    printf("M anipulation ");
                    break;
                default:
                    break;
                }
                printf("|");
                if((row+side_edge)<4) printf("                                   |");
                else printf("                                    ");
                printf("\n");
            };
            printf("\n");
            Sleep(30);
        }
    }
    system("cls");
    //This section is for the upper edge to disapper
    {
        for(int upper_edge=0;upper_edge<61;upper_edge++){
            //printf("upper edge = %d\n",upper_edge);
            printf("\033[H");
            if(upper_edge==0) printf("/*");else printf("  ");
            for(int i=0;i<upper_edge/2;i++){printf(" ");}
            if(upper_edge!=61) for(int i=0;i<(60-upper_edge);i++){printf("-");}
            for(int i=0;i<upper_edge/2;i++){printf(" ");}
            if(upper_edge==0) printf("*\\"); else printf("  ");
            printf("\n");

            //Now we will add the equal to sign for the top if the triangle
            printf(" ");
            for(int equal=0;equal<10;equal++) printf("=");
            for(int equal=0;equal<14;equal++) printf(" ");
            printf(" |");
            printf("                                    ");
            printf("\n");

            
            //There are a total of 5 lines (rows) and this will be noted with the variable row
            
            for(int row=0;row<4;row++)
            {
                //Now is the first row we will decided how many spaces must be added before printing the traingle
                for(int space=0;space<row+1;space++){printf(" ");}
                //And then we will print the double edge of the trianlge
                printf("\\\\");
                //then depending on the row, we will again add the space
                for(int space=0;space<2*(4-row)-1;space++){printf(" ");}
                //And then we print the other edge of the trinlge
                printf("/");
                //After the trianlge is printed we will again print the space
                for(int space=0;space<row+1;space++){printf(" ");}
                //Now we type the open form full form
                switch (row)
                {
                case 0:
                    printf("F ield        ");
                    //Now we may add the closing sign
                    break;
                case 1:
                    printf("O peration    ");
                    break;
                case 2:
                    printf("A nd          ");
                    break;
                case 3:
                    printf("M anipulation ");
                    break;
                default:
                    break;
                }
                printf("|");
                printf("\n");
            };
            
            //And now the letters that is FOAM will move to the sides and will start to blink

            Sleep(25);
        }
    }
    system("cls");
    //This section will move the letter to the right
    {
        for(int letter_move=0;letter_move<15;letter_move++){
            printf("\033[H");printf("\033[K");
            printf("\n");
            //Now we will add the equal to sign for the top if the triangle
            printf(" ");
            for(int equal=0;equal<10;equal++) printf("=");
            for(int equal=0;equal<14;equal++) printf(" ");
            for(int space=0;space<letter_move*5;space++){printf(" ");} printf(" |");
            //printf("                                    ");
            printf("\n");

            
            //There are a total of 5 lines (rows) and this will be noted with the variable row
            for(int row=0;row<4;row++)
            {
                //Now in the first row we will decide how many spaces must be added before printing the traingle
                for(int space=0;space<row+1;space++){printf(" ");}
                //And then we will print the double edge of the trianlge
                printf("\\\\");
                //then depending on the row, we will again add the space
                for(int space=0;space<2*(4-row)-1;space++){printf(" ");}
                //And then we print the other edge of the trinlge
                printf("/");
                //After the trianlge is printed we will again print the space
                for(int space=0;space<row+1;space++){printf(" ");}
                
                //Now we type the open form full form
                switch (row)
                {
                    case 0:
                        for(int space=0;space<letter_move;space++){printf(" ");}
                        printf("F ield        ");
                        for(int space=0;space<letter_move*4;space++) printf(" ");
                        printf("|");
                        //Now we may add the closing sign
                        break;
                    case 1:
                        for(int space=0;space<letter_move*2;space++){printf(" ");}
                        printf("O peration    ");
                        for(int space=0;space<letter_move*3;space++) printf(" ");
                        printf("|");
                        break;
                    case 2:
                        for(int space=0;space<letter_move*3;space++){printf(" ");}
                        printf("A nd          ");
                        for(int space=0;space<letter_move*2;space++) printf(" ");
                        printf("|");
                        break;
                    case 3:
                        for(int space=0;space<letter_move*4;space++){printf(" ");}
                        printf("M anipulation ");
                        for(int space=0;space<letter_move;space++) printf(" ");
                        printf("|");
                        break;
                    default:
                        break;
                    
                }
                
                printf("\n");
            };
            Sleep(100);
        }
    }
    system("cls");
    printf("\033[0m");

    //Turn the letters dim : "As if they are about to die"
    {
        printf("\033[H");
        printf("\033[2m\n");
        //Now we will add the equal to sign for the top if the triangle
        printf(" ");
        for(int equal=0;equal<10;equal++) printf("=");
        for(int equal=0;equal<14;equal++) printf(" ");
        for(int space=0;space<14*5;space++){printf(" ");} printf(" |");
        printf("                                    ");
        printf("\n");

        
        //There are a total of 5 lines (rows) and this will be noted with the variable row
        for(int row=0;row<4;row++)
        {
            //Now in the first row we will decide how many spaces must be added before printing the traingle
            for(int space=0;space<row+1;space++){printf(" ");}
            //And then we will print the double edge of the trianlge
            printf("\\\\");
            //then depending on the row, we will again add the space
            for(int space=0;space<2*(4-row)-1;space++){printf(" ");}
            //And then we print the other edge of the trinlge
            printf("/");
            //After the trianlge is printed we will again print the space
            for(int space=0;space<row+1;space++){printf(" ");}
            
            //Now we type the open form full form
            switch (row)
            {
                case 0:
                    for(int space=0;space<14;space++){printf(" ");}
                    printf("F ield        ");
                    for(int space=0;space<14*4;space++) printf(" ");
                    printf("|");
                    //Now we may add the closing sign
                    break;
                case 1:
                    for(int space=0;space<14*2;space++){printf(" ");}
                    printf("O peration    ");
                    for(int space=0;space<14*3;space++) printf(" ");
                    printf("|");
                    break;
                case 2:
                    for(int space=0;space<14*3;space++){printf(" ");}
                    printf("A nd          ");
                    for(int space=0;space<14*2;space++) printf(" ");
                    printf("|");
                    break;
                case 3:
                    for(int space=0;space<14*4;space++){printf(" ");}
                    printf("M anipulation ");
                    for(int space=0;space<14;space++) printf(" ");
                    printf("|");
                    break;
                default:
                    break;
                
            }
            printf("\n");
        };
        printf("\033[0m"); 
    }
    Sleep(50);
    system("cls");
    //This will make the FOR OUR ACADEMIC MISERY part
    {
        for(int frame = 0;frame<6;frame++)
        {   
            //printf("frame = %d\n",frame);
            printf("\033[H");
            printf("\n");

            //Now we will add the equal to sign for the top if the triangle
            
            if(frame==0)
            {
                printf(" ");
                for(int equal=0;equal<10;equal++) printf("=");
                for(int equal=0;equal<14;equal++) printf(" ");
            }
            else{
                for(int i=0;i<4;i++)printf(" ");
                printf("\033[1;32m");
                printf("//");
                printf(" \\");
                printf("\033[0m");
            }
            
            //For the first frame, the bomb will appear to be yellow and bright
            if(frame==0){printf("    \033[93mor            ur            cademic         isery\033[0m                  |");}
            else {for(int space=0;space<14*5+1;space++){printf(" ");} printf("\033[1;32m                 |\033[0m");}
            printf("\n");

            //There are a total of 5 lines (rows) and this will be noted with the variable row
            for(int row=0;row<4;row++)
            {
                //printf("row = %d\n",row);
                if(frame>1&&row==0) printf("\033[32m   //   \\   \033[0m");
                else if(row==0) printf(" \\\\       /");
                if(frame>2&&row==1) printf("\033[32m  //     \\ \033[0m");
                else if(row==1)printf("  \\\\     / ");
                if(frame>3&&row==2) printf("\033[32m //       \\\033[0m");
                else if(row==2)printf("   \\\\   / ");
                if(frame>4&&row==3) printf("\033[32m ========== \033[0m");
                else if(row==3)printf("    \\\\ /  ");

                
                if(row == 0)
                {                       
                    if(frame == 0){      printf("               F ield                                                                |");}
                    else if(frame==1){   printf("               F \033[93mor\033[0m            \033[93mur\033[0m            \033[93mcademic\033[0m         \033[93misery\033[0m                   |");}
                    else if(frame==2){   printf("               F \033[93mor\033[0m                                                                 |");}
                    else if(frame==3){   printf("               F \033[93mor\033[0m                                                                 |");}
                    else if(frame==4){   printf("               F \033[93mor\033[0m                                                                 |");}
                    else{               printf("               F \033[93mor\033[0m                                                                 |");}
                }
                else if(row==1)
                {
                    if(frame == 0){      printf("                             O peration                                              |");}
                    else if(frame==1){   printf("                 \033[2;31mield\033[0m        O peration                                              |");}
                    else if(frame==2){   printf("                             O \033[93mur\033[0m            \033[93mcademic\033[0m         \033[93misery\033[0m                   |");}
                    else if(frame==3){   printf("                             O \033[93mur\033[0m                                                    |");}
                    else if(frame==4){   printf("                             O \033[93mur\033[0m                                                    |");}
                    else{               printf("                             O \033[93mur\033[0m                                                    |");}
                }
                else if(row==2)
                {
                    if(frame == 0){      printf("                                           A nd                                       |");}
                    else if(frame==1){   printf("                                           A nd                                       |");}
                    else if(frame==2){   printf("                 \033[2;31meild\033[0m          \033[2;31mperation\033[0m    A nd                                       |");}
                    else if(frame==3){   printf("                                           A \033[93mcademic\033[0m         \033[93misery\033[0m                    |");}
                    else if(frame==4){   printf("                                           A \033[93mcademic\033[0m                                 |");}
                    else{               printf("                                           A \033[93mcademic\033[0m                                 |");}
                }
                else if(row==3)
                {
                    if(frame == 0){      printf("                                                           M anipulation              |");}
                    else if(frame==1){   printf("                                                           M anipulation              |");}
                    else if(frame==2){   printf("                                                           M anipulation              |");}
                    else if(frame==3){   printf("                 \033[2;31meild\033[0m          \033[2;31mperation\033[0m      \033[2;31mnd\033[0m            M anipulation              |");}
                    else if(frame==4){   printf("                                                           M \033[93misery\033[0m                    |");}
                    else{               printf("                                                           M \033[93misery\033[0m                  |");}
                }
                
                printf("\n");
            };
            printf(" ");
            for(int equal=0;equal<10;equal++) printf(" ");
            for(int equal=0;equal<14;equal++) printf(" ");
            if(frame==4){                printf("    \033[2;5;31mield          peration      nd              anipulation            |\033[0m");Sleep(100);}
            printf("\n");
            Sleep(5000);
        }
    }    
    
    //Now it's time to rearrange everything in the same position.
    {
        for(int frame=0;frame<=15;frame++)
        {
            //printf("frame = %d\n",frame);
            printf("\033[H");
            system("cls");
            printf("\n"); //the top space
            //And now we print the rows
            for (int row=0;row<5;row++)
            {
                if(row!=4)
                {
                    //This will make the triangle ready.
                    for(int i=0;i<4-row;i++) printf(" ");
                    printf("\033[32m//\033[0m");
                    
                    for(int i=0;i<row*2+1;i++)printf(" ");
                    printf("\033[32m\\\033[0m");
                    for(int i=0;i<4-row;i++) printf(" ");
                    if(row==0){for(int i=0;i<84-frame*5;i++){printf(" ");}printf("\033[32m|\033[0m");}
                }
                //And now we will bring the letters in the frame back

                switch (row)
                {
                    case 1:
                        for(int i=0;i<15-frame;i++) printf(" ");
                        printf("F \033[1;33mor\033[0m");
                        for(int i=0;i<65-frame*4;i++) printf(" ");printf("|");
                        break;
                    case 2:
                        for(int i=0;i<30-frame*2;i++) printf(" ");
                        printf("O \033[1;33mur\033[0m");
                        for(int i=0;i<50-frame*3;i++) printf(" ");printf("|");
                        break;
                    case 3:
                        for(int i=0;i<45-frame*3;i++) printf(" ");
                        printf("A \033[1;33mcademic\033[0m");
                        for(int i=0;i<30-frame*2;i++) printf(" ");printf("|");
                        break;
                    case 4:                 
                        //This will print the final eqaul to sign
                        printf(" ");for(int i=0;i<10;i++){printf("\033[32m=\033[0m");}
                        printf("  ");
                        //This is the lie containing misery
                        //Space before M isery (58 spaces)
                        for(int i=0;i<60-frame*4;i++) printf(" ");
                        printf("M \033[1;33misery\033[0m");
                        for(int i=0;i<16-frame*1;i++) printf(" ");printf("|");
                        break;
                    default:
                        break;
                }
                printf("\n");
            }
            Sleep(50);
        }
    }
    
    printf("\033[?25h");//Show the curser again
    return;
}

void logo(){
    return;
}
/* INPUT AND OUTPUT ********************************************************************************/

/*Now is the time we will have to ask for the type of problem to be solve. At first the options are
are limirted to heat equation but the later studies will include involving other model as well*/
void input_type_of_problem()
{
    printf("\033[1;33m");
    printf("Which type of problem do you want to solve: \n");
    printf("1. Heat Equation \n"
        "2. Laplace Equation \n"
        "3. Poisson Equation\n"
        "4. Wave Equation\n"
        "5. Advection-Diffusion Equation\n"
        "6. Burger Equation\n"
        "7. Euler-Bernoulli Beam Equation\n"
        "8. Diffusion Equation\n"
        "9. Elasticity Equations\n"
        "10. Reaction-Diffusion Equations\n"
        "11. Population dynamics PDEs\n"
        "12. Navier-Stokes Equations\n"
        "13. Radiation Heat Transfer Equation\n"
        "14. Laplace-Young Equation\n");
        printf("\nEnter the problem number: ");
        scanf("%d",&problem_type);
    system("cls");
    printf("Whats the dimension of problem: ");
    scanf("%d",&dimension);
    printf("What kind of scheme (Explicit: 1| Implicit: 2): ");
    scanf("%d",&scheme_type);
    printf("\033[0m\n");
    
    if(scheme_type==1) // If the scheme is explicit
    {
        if(problem_type==1) // Heat Eqaution
        {
            printf("You have selected \033[1m Explicit scheme\033[0m for \033[1m Heat equation\033[0m\n");
            printf("Here are the list of schemes to choose from\n"
                "1. Forward in Time and Central in Space (FTCS)\n"
                "2. Dufort Frankel Scheme\n");
            printf("Enter the scheme to choose: ");
            scanf("%d",&scheme);
        }
        else if(problem_type==2) // Laplace Equation
        {
            printf("You have selected \033[1m Explicit scheme\033[0m for \033[1m Laplace equation\033[0m\n");
            printf("Here are the list of schemes to choose from\n"
                "1. Successive Over-Relaxation (SOR) Method\n");
            printf("Enter the scheme to choose: ");
            scanf("%d",&scheme);
        }
        else if(problem_type==4) // Wave Equation
        {
            printf("You have selected \033[1m Explicit scheme\033[0m for \033[1m Wave equation\033[0m\n");
            printf("\033[4;33mLinear Wave Equation\033[0m\n");
            printf("Here are the list of schemes to choose from\n"
                "1. Euler FTFS\n"
                "2. Euler FTCS\n"
                "3. First order upwiding scheme\n"
                "4. Lax Method\n"
                "5. Lax Wendroiff Scheme\n"
                "6. Lax Wendroff MultiStep method\n"
                "7. Mac Cormack method\n"
            );
            printf("\n\033[4;33mNon-Linear Wave Equation\033[0m\n");
            printf("Here are the list of schemes to choose from\n"
                    "1. Lax Method\n"
                    "2. Lax Wendroff Method\n"
                    "3. MacCormack Method\n"
            );
            printf("Enter the scheme to choose: "); scanf("%d",&scheme);
        }
    }

    else if(scheme_type==2) //If the scheme is implicit
    {
        if(problem_type==1) //Heat Eqaution
        {
            printf("You have selected \033[1m Implicit scheme\033[0m for \033[1m Heat equation\033[0m\n");
            printf("Here are the list of schemes to choose from\n"
                "1. Backward Euler (Implicit Euler) Method\n"
                "2. Crank-Nicolson Method\n");
            printf("Enter the scheme to choose: ");
            scanf("%d",&scheme);
        }
        else if(problem_type==2) //Laplace Equation
        {
            printf("You have selected \033[1m Implicit scheme\033[0m for \033[1m Laplace equation\033[0m\n");
            printf("Here are the list of schemes to choose from\n"
                "1. SOR by line\n");
            printf("Enter the scheme to choose: ");
            scanf("%d",&scheme);
        }
        else if(problem_type==4) //Wave Eqaution
        {
            printf("You have selected \033[1m Implicit scheme\033[0m for \033[1m Wave equation\033[0m\n");
            printf("\033[4;33mLinear Wave Equation\033[0m\n");
            printf("Here are the list of schemes to choose from\n"
                "1. Euler's BTCS Method\n"
                "2. Crank Nicholsen Method\n"
            );
            printf("\033[4;33mNon-Linear Wave Equation\033[0m\n");
            printf("Here are the list of schemes to choose from\n"
                "1. Beam Warming Implicit Method\n"
            );
            printf("Enter the scheme to choose: ");
            scanf("%d",&scheme);
        }
    }
    return;
}

/*This function will ask for the input finally. And as decided the input will be first to ask for the 
first matrix which will be the values, and then the next will be to ask for the boudary matrix*/
input_1 input_data_1(input_1 input){
    printf("What is the length of domain: "); scanf("%lf",&L_x);
    //Then we will ask for the grid points
    printf("Enter the grid points: ");
    scanf("%d",&grid_points_x);
    
    //The first index will tell which matrix is it? Value (0) or Boundary(1)
    input.boundary_matrix = (int*)malloc(sizeof(int)*grid_points_x);
    input.value_matrix = (double*)malloc(sizeof(double)*grid_points_x);
    
    //Now is the time to scan the input
    printf("Enter the Values Matirx (First Matrix)\n");
    for(int i=0;i<grid_points_x;i++){scanf("%lf",&input.value_matrix[i]);};
    printf("Values stored inside the values matrix are: \n");
    for(int i=0;i<grid_points_x;i++){printf("%lf ",input.value_matrix[i]);};

    printf("\nInstruction to enter Boundary Matrix:\n"
        "1: Forced Dircihlet Boundary Condition 1\n"
        "2: Forced Dirichlet Boundary Condition 2\n"
        "3: Forced Dirichlet Boundary Conditon 3\n"
        "4: Forced Neumann Boundary Condition 1\n"
        "5: Forced Neumann Boundary Condition 2\n"
        "6: Forced Neumann Boundary Condition 3\n"
        "7: Custom Boundary Condition\n");
    printf("Enter the Boudary Matrix (Second Matrix)\n");
    for(int i=0;i<grid_points_x;i++){scanf("%d",&input.boundary_matrix[i]);}; 
    printf("Values stored inside the boundary matrix are: \n");
    for(int i=0;i<grid_points_x;i++){printf("%d ",input.boundary_matrix[i]);};
    printf("\n");
    return input;
}

input_2 input_data_2(input_2 input)
{
    //We'll ask for lenght of domain
    printf("What is the length of domain in X: "); scanf("%lf",&L_x);
    printf("What is the length of domain in Y: "); scanf("%lf",&L_y);
    //Then we will ask for the grid points
    printf("Enter the grid points X: "); scanf("%d",&grid_points_x);
    printf("Enter the grid points Y: "); scanf("%d",&grid_points_y);

    //This means the first index will be x cordinate
    input.boundary_matrix = (int**)malloc(sizeof(int*)*grid_points_x);
    input.value_matrix = (double**)malloc(sizeof(double)*grid_points_x);
    //This means the second index is that of y cordinate
    for(int k=0;k<grid_points_x;k++){input.boundary_matrix[k]= (int*)malloc(sizeof(int)*grid_points_x);input.value_matrix[k] = (double*)malloc(sizeof(double)*grid_points_x);};
    for(int k=0;k<grid_points_y;k++){input.boundary_matrix[k] = (int*)malloc(sizeof(int)*grid_points_y);input.value_matrix[k] = (double*)malloc(sizeof(double)*grid_points_y);};
    
    //Now is the time to scan the input
    /*Now here is something important. I am considering the point (x,y) = (0,0) to be located at the bottom left
    and the point (L_x,L_y) to be at the top right*/
    printf("Enter the Values Matirx (First Matrix)\n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            scanf("%lf",&input.value_matrix[i][j]);
        }
    };
    printf("Values stored inside the values matrix are: \n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            printf("%lf ",input.value_matrix[i][j]);
        }; printf("\n");};

    printf("Instruction to enter Boundary Matrix:\n"
        "1: Forced Dircihlet Boundary Condition 1\n"
        "2: Forced Dirichlet Boundary Condition 2\n"
        "3: Forced Dirichlet Boundary Conditon 3\n"
        "4: Forced Neumann Boundary Condition 1\n"
        "5: Forced Neumann Boundary Condition 2\n"
        "6: Forced Neumann Boundary Condition 3\n"
        "7: Custom Boundary Condition\n");
    printf("Enter the Boundary Matrix (Second Matrix)\n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            scanf("%d",&input.boundary_matrix[i][j]);
        }
    };
    printf("Values stored inside the boundary matrix (Second Matrix)\n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            printf("%d ",input.boundary_matrix[i][j]);
        }printf("\n");
    };
    return input;
}

void output_1(input_1 input)
{
    printf("Into the function output1D\n");
    if(problem_type==1)//This means the problem is Heat Equation
    {
        printf("So you have selected to solve the heat equation\n");
        Heat_Eqaution_1(input);
        printf("The program is executed\n");
    }
    else if(problem_type==2)//This means the problem is Laplace Equation
    {
        printf("So you have selected to solve the Laplace equation\n");
        Laplace_Eqaution_1(input);
        printf("The program is executed\n");
    }
    else if(problem_type==4) //This means the problem is Wave Equation
    {
        int islinear;
        printf("Linear Wave Eqaution (1) OR Non-Linear Wave Eqaution/Invicid Burger's eqaution (0) : ");
        scanf("%d",&islinear);
        printf("So you have selected to solve the Wave equation\n");
        input = Wave_Eqaution_1(input,islinear);
        printf("The program is executed\n");
    }
    print_output_1(input);
    return;
}

void output_2(input_2 input)
{
    printf("Into the function output 2D\n");
    if(problem_type==1)//This means the problem is Heat Equation
    {
        printf("So you have selected to solve the Heat equation\n");
        Heat_Eqaution_2(input);
        
    }
    else if(problem_type==2)//This means the problem is Laplace Equation
    {
        printf("So you have selected to solve the Laplace equation\n");
        Laplace_Eqaution_2(input);
        
    }
    print_output_2(input);
    return;
}

/* STREAM FUNCTION AND VORTICITY FORMATIONS **********************************************************/

//These are some important global variable
double dynamic_viscocity;
double characteristic_lenght;
double characteristic_velocity;

//Input contains two matrix, first one being the X-component matrix, and the second matrix being Y-velocity matrix.
typedef struct input_data_stream_function_vorticity {
    double ** u;// This matrix will be to store x velocity
    double ** v;// This matrix will be to store y velocity
    int ** boundary_matrix; //This will be for the boundary condition
    double ** w; //This will store the vorticity at each point
    double ** s;// This will stroe the stream function at each point
    double ** p; //This will store the pressure at each point
}input_stream_vorticity_2;

//This will set up the input matrix from the user (Custom input)
input_stream_vorticity_2 input_data_stream_function_vorticity_2(input_stream_vorticity_2 input)
{
    printf("Enter lenght of domain in X direction: "); scanf("%d",&L_x);
    printf("Enter lenght of domain in Y direction: "); scanf("%d",&L_y);
    printf("Enetr grid size in X direction: "); scanf("%d",&grid_points_x);
    printf("Enetr grid size in Y direction: "); scanf("%d",&grid_points_y);

    //Now we allocate the memory
    //First index will be x cordinate
    // Second index will be y cordinate
    input.boundary_matrix = (int**)malloc(sizeof(int*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.boundary_matrix[i] = (int*)calloc(sizeof(int),grid_points_y);}
    input.u = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.u[i] = (double*)calloc(sizeof(double),grid_points_y);}
    input.v = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.v[i] = (double*)calloc(sizeof(double),grid_points_y);}
    input.p = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.p[i] = (double*)calloc(sizeof(double),grid_points_y);}
    input.w = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.w[i] = (double*)calloc(sizeof(double),grid_points_y);}
    input.s = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.s[i] = (double*)calloc(sizeof(double),grid_points_y);}
    //So it is imporat to mention that the values are stored in the cartesian cordinate based way. It starts for x from 0 to X, and y from Y to 0, as the input is typed.
    //Now we take the input from the user
    printf("Enter the X Velocities Matrix\n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            scanf("%lf",&input.u[i][j]);
        }
    };
    printf("Enter the Y Velocities Matrix\n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            scanf("%lf",&input.v[i][j]);
        }
    };
    //Now we print the values see validate it
    printf("Values stored inside the X velocity matrix are: \n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            printf("%lf ",input.u[i][j]);
        }; printf("\n");};
    //Now we print the values see validate it
    printf("Values stored inside the Y velocity matrix are: \n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            printf("%lf ",input.v[i][j]);
        }; printf("\n");};

    //So it is interesting to make the boundary condition work. We first have the following:
    /*So here are the revised instruction for setting the boundary points in the domain.The domian will always be a rectangle, and point (0,0)
    will be assigned the value of stream function equal to zero. And then all the values of stream fucntion along the boundary will be calculated
    accordign to this point's reference.
    Now the boundary point can have */
    printf("Instruction to Enter Boundary Matrix:\n"
        "1: Forced Dircihlet Boundary Condition 1\n"
        "2: Forced Dirichlet Boundary Condition 2\n"
        "3: Forced Dirichlet Boundary Conditon 3\n"
        "4: Forced Neumann Boundary Condition 1\n"
        "5: Forced Neumann Boundary Condition 2\n"
        "6: Forced Neumann Boundary Condition 3\n"
        "7: Custom Boundary Condition\n");
    printf("Enter the Boundary Matrix\n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            scanf("%d",&input.boundary_matrix[i][j]);
        }
    };
    printf("Values stored inside the boundary matrix\n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++){
            printf("%d ",input.boundary_matrix[i][j]);
        }printf("\n");
    };

    //So now the input matrix is filled with X velocities and Y velociites and Boundary values.
    return input;
}

// A simple lid cavity problem. (Test input)
input_stream_vorticity_2 test_input_function(input_stream_vorticity_2 input,const char** version,const char** input_type,int choice_version,int input_type_number)
{
    system("cls");
    logo();
    printf("\033[1;34m\nVersion %d\t:\t %s\033[0m\n",choice_version+1,version[choice_version]);
    printf("\033[1;34mProblem Type\t:\t %s\033[0m\n\n",input_type[input_type_number]);

    printf("Enter lenght of domain in X direction: "); scanf("%lf",&L_x);
    printf("Enter lenght of domain in Y direction: "); scanf("%lf",&L_y);
    printf("Enter grid size in X direction: "); scanf("%d",&grid_points_x);
    printf("Enter grid size in Y direction: "); scanf("%d",&grid_points_y);

    //Now we allocate the memory
    //First index will be x cordinate
    // Second index will be y cordinate
    input.boundary_matrix = (int**)malloc(sizeof(int*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.boundary_matrix[i] = (int*)calloc(sizeof(int),grid_points_y);}
    input.u = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.u[i] = (double*)calloc(sizeof(double),grid_points_y);}
    input.v = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.v[i] = (double*)calloc(sizeof(double),grid_points_y);}
    input.p = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.p[i] = (double*)calloc(sizeof(double),grid_points_y);}
    input.w = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.w[i] = (double*)calloc(sizeof(double),grid_points_y);}
    input.s = (double**)malloc(sizeof(double*)*grid_points_x); for(int i=0;i<grid_points_x;i++){input.s[i] = (double*)calloc(sizeof(double),grid_points_y);}
    
    if(debugger) printf("\033[32mMemory is allocated\033[0m\n");
    //Here the x velocity is assigned the value 10 m/s in the right direction (positive), and y velocity to be zero.
    for(int i=0;i<grid_points_x;i++)
    {
        if(i>0 && i<grid_points_x-1) input.u[i][grid_points_y-1] = 10;
        input.v[i][grid_points_y-1] = 0;
    }
    if(debugger) printf("\033[32mValues at top the boundary\033[0m\n");
    //And rest of all the points are assigned the velocity of 0 m/s
    for(int i=0;i<grid_points_x;i++)
        {
            for(int j=0;j<grid_points_y-1;j++)
            {
                input.u[i][j] = 0; input.v[i][j] = 0;
            }
        }
    if(debugger) printf("\033[32mValues of Y velocity is allocated\033[0m\n");
    //For now there is no need of boundary since the program doesnt handles the boundary values for now.
    //So now the input matrix is filled with X velocities and Y velociites and Boundary values.
    //Here is the domain and the values of u and v assigned to them
    if(debugger) //Prints the velocity matrices
    {    
        printf("\033[1m");
        printf("X velocity Matrix : \n");
        printf("\033[32m");
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0 ;i<grid_points_x;i++)
            {
                printf("%lf ",input.u[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
        printf("\033[1m Y Velocity Matrix : \n");
        printf("\033[32m");
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0 ;i<grid_points_x;i++)
            {
                printf("%lf ",input.v[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
    }
    
    char proceed;
    printf("You want to proceed (yes (y) | no (n))? "); 
    scanf(" %c",&proceed);
    //printf("So the character choosen is %c\n",proceed);

    if(proceed == 'y'){return input;}
    else {terminate = true; return input;}
   
}

//Step 0: We will first non dimensionalize the velocities:
input_stream_vorticity_2 non_dimensionalise_velocity(input_stream_vorticity_2 input,double U,double* delta_x,double* delta_y,double* delta_t,double L){
    //This function will take in the input matrix, and the chracteristic velocity to non dimensionalise the velocity
    for(int i=0;i<grid_points_x;i++){for(int j=0;j<grid_points_y;j++){
        input.u[i][j] = input.u[i][j]/U;
        input.v[i][j] = input.v[i][j]/U;
    }}
    //This will non dimensionalise the lenghts
    *delta_x = *delta_x/L;
    *delta_y = *delta_y/L;
    *delta_t = *delta_t*(characteristic_velocity/characteristic_lenght);
    return input;
    /*CONFUSION
    /After the velocities are non dimensioanlised we could move on with our calculation. Now one confusion that I get in what abotu non dimensionalising the lenghts. But I think that
    when formulating the navier strokes equation we could use the reynolds number, and that will incorporate the data about lengths and dynamics viscocity. So
    to strat with, the only thing we need is non dimensionalised velocity*/
}

//Step 1: This function will take in the input and calculate the vorticity at the interior points at time t. (Taylor Series)
input_stream_vorticity_2 velocity_at_t1_to_interior_vorticity_at_t1(input_stream_vorticity_2 input,double delta_x,double delta_y)
{
    //Welcome debugger
    if(debugger)
    {
        printf("\033[33mInto the function velocity_at_t1_to_vorticity_at_t1\n");
        printf("delta_x = %lf\n",delta_x);
        printf("delta_y = %lf\n",delta_y);
    } 
    
    //This will start the clock
    time_t start = clock();
    //Computing the values of interior vorticity using complete velocity through out the domain. 
    for(int j=grid_points_y-2;j>=1;j--)
    {
        for(int i=1;i<grid_points_x-1;i++)
        {
            //Using the formula interior vorticity is equal to 
            input.w[i][j] = (input.v[i+1][j] - input.v[i-1][j])/(2*delta_x) - (input.u[i][j+1] - input.u[i][j-1])/(delta_y);
            if(debugger && strong_debugger) //Printing all the components
            {
                printf("\033[31mi = %d | j = %d\033[0m\n",i,j);
                printf("\033[32m");
                printf("v(i+1,j) = %lf\n",input.v[i+1][j]);
                printf("v(i-1,j) = %lf\n",input.v[i-1][j]);
                printf("\n");
                printf("u(i,j+1) = %lf\n",input.u[i][j+1]);
                printf("u(i,j-1) = %lf\n",input.u[i][j-1]);
                printf("\n");
                printf("w(i,j) = %lf\n",input.w[i][j]);
                printf("\n");
                printf("\033[0m");
            }
        }
        if(debugger && strong_debugger) printf("\n");
    }
    
    //This will end the clock
    time_t end = clock();
    if(computation_time) //This will print the runtime
    {
        double cpu_time = (double)(end-start)/CLOCKS_PER_SEC;
        printf("\033[1;34m");
        printf("Interior Vorticity Calculation CPU Time : %lf secs",cpu_time);
        printf("\033[0m");
    }
    //Printing the values of vorticity
    if(debugger)
    {
        printf("\033[1m\nVorticty matrix\033[0m\n");
        printf("\033[32m");
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0;i<grid_points_x;i++)
            {
                printf("%lf ",input.w[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
    }
    return input;
}

//Step 2: Stream Function at the boundary. This has to be done manually
input_stream_vorticity_2 velocity_to_stream_function_at_boundary_at_t1(input_stream_vorticity_2 input,double delta_x,double delta_y)
{
    if(debugger){printf("\033[31mInto stream_function_at_boundary_at_t1\033[0m\n");}
    time_t start = clock();
    //This function will update the stream functions at the boundary, and this will be automated. 
    /*So here is what i know about the boundary. They are assigned to me such that */
    /*So we will start from bottom side of the domain, and to the right side, then to top and then to left side of the domain, eventually reaching the same point*/
    if(debugger)
    {
        printf("delta_x = %lf\n",delta_x);
        printf("delta_y = %lf\n",delta_y);
        printf("grid_points_x = %d\n",grid_points_x);
        printf("grid_points_y = %d\n",grid_points_y);
    }
    
    input.s[0][0] = 0;
    printf("\n");
    //Bottom Side
    for(int i=0;i<grid_points_x-1;i++)
    { 
        if(debugger && strong_debugger)
        {
            printf("v(%d,0) = %lf\n",i,input.v[i][0]);
            printf("u(%d,0) = %lf\n",i,input.u[i][0]);
            printf("s(%d,0) = %lf\n",i,input.s[i][0]);
        }
        
        //This comes after solving the equation : y-velocity equals negative of partial derivative of stream function wrt x
        if(i==0) //Corner Problem
        {
            if(input.v[i][0] && !input.v[i+1][0]) input.s[i+1][0] = input.s[i][0];
            else {input.s[i+1][0] = input.s[i][0] - input.v[i][0]*(delta_x);}
        }
        else {input.s[i+1][0] = input.s[i][0] - input.v[i][0]*(delta_x);}
        //This comes after solving the equation : x-velocity equals positive of partial derivative of stream function wrt y
        input.s[i][1] = input.s[i][0] + input.u[i][0]*(delta_y);
        
        if(debugger && strong_debugger)
        {
            printf("s(%d,1) = %lf\n",i,input.s[i][1]);
            printf("s(%d,0) = %lf\n",i+1,input.s[i+1][0]);
            printf("\n");
        }
    }
    //This will print the values of streamfucntion at the bottom side of domain
    if(debugger && strong_debugger)
    {
        printf("\033[31m");
        printf("Stream fucntion at the Bottom side of domian\n");
        for(int i=0;i<grid_points_x;i++){printf("%lf ",input.s[i][0]);}
        printf("\033[0m\n");
        printf("\n");
    }

    //Rigth side
    const int I = grid_points_x-1;
    for(int j=0;j<grid_points_y-1;j++)
    {
        if(debugger && strong_debugger)
        {
            printf("v(%d,%d) = %lf\n",I,j,input.v[I][j]);
            printf("u(%d,%d) = %lf\n",I,j+1,input.u[I][j]);
            printf("s(%d,%d) = %lf\n",I,j,input.s[I][j]);
        }
        
        //This comes after solving the equation : y-velocity equals negative of partial derivative of stream function wrt x
        if(j==0) //Corner problem
        {
            if(input.u[I][j] && input.u[I][j+1]) input.s[I][j+1] = input.s[I][j];
            else {input.s[I][j+1] = input.s[I][j] + input.u[I][j]*delta_y;}
        }
        else input.s[I][j+1] = input.s[I][j] + input.u[I][j]*delta_y;
        //This comes after solving the equation : x-velocity equals positive of partial derivative of stream function wrt y 
        input.s[I-1][j] = input.s[I][j] + input.v[I][j]*delta_x;
        if(debugger && strong_debugger)
        {
            printf("s(%d,%d) = %lf\n",I,j+1,input.s[I][j+1]);
            printf("s(%d,%d) = %lf\n",I-1,j,input.s[I-1][j]);
            printf("\n");
        }
    }
    //This will print the values of streamfucntion at the right side of domain
    if(debugger && strong_debugger)
    {
        printf("\033[31m");
        printf("Stream fucntion at the Right side of domian\n");
        for(int j=grid_points_y-1;j>=0;j--){printf("%lf\n",input.s[I][j]);}
        printf("\033[0m\n");
        printf("\n");
    }

    //Top Side
    const int J = grid_points_y-1;
    for(int i=grid_points_x-1;i>0;i--)
    {
        if(debugger && strong_debugger)
        {
            printf("u(%d,%d) = %lf\n",i,J,input.u[i][J]);
            printf("v(%d,%d) = %lf\n",i,J,input.v[i][J]);
            printf("s(%d,%d) = %lf\n",i,J,input.s[i][J]);
        }
        //This comes after solving the equation : y-velocity equals negative of partial derivative of stream function wrt x
        if(i==I) //Corner problem
        {
            if(input.v[i][J] && !input.v[I-1][J]) input.s[i-1][J] = input.s[i][J];
            else input.s[i-1][J] = input.s[i][J] + input.v[i][J] * delta_x;
        }
        else input.s[i-1][J] = input.s[i][J] + input.v[i][J] * delta_x;
        //This comes after solving the equation : x-velocity equals positive of partial derivative of stream function wrt y 
        input.s[i][J-1] = input.s[i][J] - input.u[i][J]*delta_y;
        if(debugger && strong_debugger)
        {
            printf("s(%d,%d) = %lf\n",i,J-1,input.s[i][J-1]);
            printf("s(%d,%d) = %lf\n",i-1,J,input.s[i-1][J]);
            printf("\n");
        }
    }
    //Printing the streamfunction of top side
    if(debugger && strong_debugger)
    {
        printf("\033[31m");
        printf("Stream function at the Top side of domain\n");
        for(int i=0;i<grid_points_x;i++){printf("%lf ",input.s[i][J]);}
        printf("\033[0m\n");
        printf("\n");
    }

    //Left Side
    for(int j=grid_points_y-1;j>0;j--)
    {
        if(debugger && strong_debugger)
        {
            printf("u(%d,%d) = %lf\n",0,j,input.u[0][j]);
            printf("v(%d,%d) = %lf\n",0,j,input.v[0][j]);
            printf("s(%d,%d) = %lf\n",0,j,input.s[0][j]);
        }
        //This comes after solving the equation : y-velocity equals negative of partial derivative of stream function wrt x
        input.s[1][j] = input.s[0][j] + input.v[0][j]*delta_x;
        //In all other cases 
        //only in the case of corner (0,J-1)
        if(j == J)
        {
            //if(strong_debugger) printf("Entered the corner function for top side\n");
            //Analysis this by itself,and you will know.
            if(input.u[0][j] && !input.u[0][j-1]) input.s[0][j-1] = input.s[0][j];
            else input.s[0][j-1] = input.s[0][j] - input.u[0][j]*delta_y;
        }
        else input.s[0][j-1] = input.s[0][j] - input.u[0][j]*delta_y;
        if(debugger && strong_debugger)
        {
            printf("s(%d,%d) = %lf\n",0,j-1,input.s[0][j-1]);
            printf("s(%d,%d) = %lf\n",1,j,input.s[1][j]);
            printf("\n");
        }
    }
    //Printing the streamfunction of left side
    if(debugger && strong_debugger)
    {
        printf("\033[31m");
        printf("Stream function at the left side of domain\n");
        for(int j=grid_points_y-1;j>=0;j--){printf("%lf \n",input.s[0][j]);}
        printf("\033[0m\n");
        printf("\n");
    }

    if(debugger)
    {
        printf("\033[31mStream function matrix\033[0m\n");
        printf("\033[32m");
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0;i<grid_points_x;i++)
            {
                printf("%lf ",input.s[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
    }
    time_t end = clock();
    if(computation_time)
    {
        double cpu_time = (double)(end-start)/CLOCKS_PER_SEC;
        printf("\033[34mStream Function Calculation Total Time : %lf secs\033[0m\n",cpu_time);
    }

    return input;
}

/*Step 3: This function will take in the input, and will update all the streamfunction points at time t in the interior using the interior 
vorticity at time t, and also make stream function equal to constant (0) at the body surface which is it's boundary condition.
At this point all the stream function values will be known at t1*/
//Scheme 1
input_stream_vorticity_2 successive_over_relaxation_stream_function_equation(input_stream_vorticity_2 input,double delta_x,double delta_y)
{
    double w = 0.5;
    double iterations = 50;
    double beta = delta_x/delta_y;

    //Start the clock
    clock_t start = clock();
    //Can also get the time taken for each iterations
    double* iteration_time = (double*)calloc(iterations,sizeof(double));
    clock_t i_start;
    clock_t i_end;
    
    for(int iteration=0;iteration<iterations;iteration++)
    {
        if(debugger) {printf("\033[31mIteration = %d\033[0m\n",iteration);}
        if(computation_time) i_start = clock();

        //Now the code will traverse from the bottom left corner to the top right corner. And here it is assumed that the boundary points are given to us. So all the points at the boundary is given
        for(int j=2;j<grid_points_y-2;j++)
        {
            for(int i=2;i<grid_points_x-2;i++)
            {
                //If the interior point doesnt have any dirichlet condition then you may update it's value
                if(!input.boundary_matrix[i][j])
                {
                    input.s[i][j] = (1/(2*(1+pow(beta,2))))*(pow(delta_x,2)*input.w[i][j] + input.s[i+1][j] + input.s[i-1][j] + pow(beta,2)*(input.s[i][j+1] + input.s[i][j-1]));
                }
            } 
        }
        if(computation_time) 
        {
            i_end = clock();
            iteration_time[iteration] = (double)(i_end-i_start)/CLOCKS_PER_SEC;
        }

        //This debugger will be used to print the matrix at every iteration.
        if(debugger && strong_debugger)
        {
            printf("\033[31mCompletion of iteration %d\033[0m\n",iteration);
            printf("\033[32m");
            for(int j=grid_points_y-1;j>=0;j--){
                for(int i=0;i<grid_points_x;i++)
                {
                    printf("%lf ",input.s[i][j]);
                }; 
                printf("\n");
            };
            printf("\033[0m\n");
        }    
    }
    //End the Clock
    clock_t end = clock();
    if(computation_time)
    {
        long double cpu_time = ((long double)(end-start)/CLOCKS_PER_SEC);
        //For iteration time
        for(int iteration = 0;iteration<iterations;iteration++)
        {
            printf("Time taken for iteration %d : %lf secs\n",iteration+1,iteration_time[iteration]);
        }
        printf("Total Time Taken : \n",cpu_time);
        free(iteration_time);
    }
    return input;
}
//Scheme 2
input_stream_vorticity_2 Thomas_Algorithm_SOR_by_line_Stream_Function_i(input_stream_vorticity_2 input,int k,double w,double beta)
{
    int n = grid_points_x-2;
    double** A = (double**)malloc(sizeof(double*)*(n));
    for(int i=0;i<n;i++){A[i] = (double*)calloc(n,sizeof(double));}
    for(int row=0;row<n;row++){
        A[row][row] = 1;
        if(row>0) A[row][row-1] = -(w/(2*(1+pow(beta,2))));
        else if(row<n-1) A[row][row+1] = -(w/(2*(1+pow(beta,2))));
    }

    double* b = (double*)malloc(sizeof(double)*n);
    /*Changed this*/
    for(int row=0;row<n;row++){b[row] = ((pow(beta,2)*w)/(2*(1+pow(beta,2))))*(input.s[k][row+2]+input.s[k][row]+input.w[k][row+1])+(1-w)*(input.s[k][row+1]);}
    
    double temp;
    for(int row=0;row<n-1;row++)
    {
        //First we make the diagonal element equal to 1
        temp= A[row][row];
        for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
        b[row] = b[row]/temp;
        
        // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
        temp = A[row+1][row];
        for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
        b[row+1] = b[row+1] - temp*b[row];
        
    }
    
    input.s[k][n] = b[n-1]/A[n-1][n-1];
    for(int row=n-1;row>=1;row--)
    {
        input.s[k][row] = b[row-1]-A[row-1][row]*input.s[k][row+1];
    }
    free(A);
    free(b);
    return input;
}
input_stream_vorticity_2 Thomas_Algorithm_SOR_by_line_Stream_Function_j(input_stream_vorticity_2 input,int k,double w,double beta)
{
    int n = grid_points_y-2;
    double** A = (double**)malloc(sizeof(double*)*(n));
    for(int i=0;i<n;i++){A[i] = (double*)calloc(n,sizeof(double));}
    
    for(int row=0;row<n;row++)
    {
        A[row][row] = 1;
        if(row>0) A[row][row-1] = -((w*pow(beta,2))/(2*(1+pow(beta,2))));
        else if(row<n-1) A[row][row+1] = -((w*pow(beta,2))/(2*(1+pow(beta,2))));
    }

    double* b = (double*)malloc(sizeof(double)*n);
    for(int row=0;row<n;row++){b[row] = (w/(2*(1+pow(beta,2))))*(input.s[row+2][k]+input.s[row][k]+input.w[row+1][k])+(1-w)*(input.s[row+1][k]);}
    
    double temp;
    for(int row=0;row<n-1;row++)
    {
        //First we make the diagonal element equal to 1
        temp= A[row][row];
        for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
        b[row] = b[row]/temp; 
        // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
        temp = A[row+1][row];
        for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
        b[row+1] = b[row+1] - temp*b[row];   
    }
    
    //Backward substitution
    input.s[n][k] = b[n-1]/A[n-1][n-1];
    for(int row=n-1;row>=1;row--)
    {
        input.s[row][k] = b[row-1]-A[row-1][row]*input.s[row+1][k];
    }
    
    free(A);
    free(b);
    return input;
}
input_stream_vorticity_2 SOR_by_line_stream_function_equation(input_stream_vorticity_2 input,double delta_x,double delta_y)
{
    if(debugger) printf("Into SOR by line function\n");
    double beta = delta_x/delta_y;
    double relaxation_factor = 0.5;
    double iterations = 50;
    for(int d = 0;d<iterations;d++)
    {
        for(int i = 1;i<grid_points_x-1;i++)
        {
            input = Thomas_Algorithm_SOR_by_line_Stream_Function_i(input,i,relaxation_factor,beta);
        }
        for(int j = 1;j<grid_points_y-1;j++)
        {
            input = Thomas_Algorithm_SOR_by_line_Stream_Function_j(input,j,relaxation_factor,beta);  
        }    
    }
    return input;
}

//Stream Function Eqaution Solver
input_stream_vorticity_2 interior_vorticity_at_t1_and_boundary_stream_at_t1_to_streamfunction_at_t1(input_stream_vorticity_2 input,double delta_x,double delta_y,int scheme)
{
    
    if(debugger){
        printf("Into Step 3\n");
        printf("Scheme = %d\n",scheme);
    }
    //Stream function eqaution (Poisson's equation)
    switch (scheme)
    {
        case 1:
            /*Selected Succesive Over Relaxation - available*/
            input = successive_over_relaxation_stream_function_equation(input,delta_x,delta_y);
            break;
        case 2:
            /*Selected SOR by line - available*/
            input = SOR_by_line_stream_function_equation(input,delta_x,delta_y);
            break;
        default:
            break;
    }
    //This will print the final stream function
    if(debugger)
    {
        printf("\033[31mFinal Stream Function\033[0m\n");
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0;i<grid_points_x;i++)
            {
                printf("%lf ",input.s[i][j]);
            }
            printf("\n");
        }
    }

    return input;
}

/*Step 4: Recompute the velocites and update them to ensure divergence free solution using the stream function ast t1*/
input_stream_vorticity_2 stream_function_at_t1_to_velocity_at_t1(input_stream_vorticity_2 input,double delta_x, double delta_y)
{
    //Velocity at the boundary is computed with second order difference scheme
    /*
        //Bottom side
        for(int i=0;i<grid_points_x;i++)
        {
            //X velocity will be calculate with forward difference and Y velocity will be calculated using central diffrence
            if(i == 0)
            {
                //This will use the second order forward difference scheme
                input.v[i][0] = -(-3*input.s[i][0] + 4*input.s[i+1][0] - input.s[i+2][0])/(2*delta_x);
            }
            else if(i == grid_points_x-1)
            {
                //This will use the second order backward difference scheme
                input.v[i][0] = -(3*input.s[i][0] - 4*input.s[i-1][0] + input.s[i-2][0])/(2*delta_x);
            }
            else 
            {
                input.v[i][0] = -(input.s[i+1][0] - input.s[i-1][0])/(2*delta_x);
            }
            input.u[i][0] = (-3*input.u[i][0] + 4*input.u[i][1] - input.u[i][2])/(2*delta_y);
        }

        //Right side
        for(int j=0;j<grid_points_y;j++)
        {
            //X velocity will be calculate with forward difference and Y velocity will be calculated using central diffrence
            if(j == 0)
            {
                //This will use the second order forward difference scheme
                input.u[grid_points_x-1][j] = (-3*input.s[grid_points_x-1][j] + 4*input.s[grid_points_x-1][j+1] - input.s[grid_points_x-1][j+2])/(2*delta_y);
            }
            else if(j == grid_points_y-1)
            {
                //This will use the second order backward difference scheme
                input.u[grid_points_x-1][j] = (3*input.s[grid_points_x-1][j] - 4*input.s[grid_points_x-1][j-1] + input.s[grid_points_x-1][j-2])/(2*delta_y);
            }
            else 
            {
                input.u[grid_points_x-1][j] = (input.s[grid_points_x-1][j+1] - input.s[grid_points_x-1][j-1])/(2*delta_y);
            }
            input.v[grid_points_x-1][j] = - (-3*input.u[grid_points_x-2][j] + 4*input.u[grid_points_x-1][j] - input.u[grid_points_x-1][j])/(2*delta_x);
        }
        
        //Top Side
        for(int i=grid_points_x-1;i>=0;i--)
        {
            //X velocity will be calculate with forward difference and Y velocity will be calculated using central diffrence
            if(i == 0)
            {
                //This will use the second order forward difference scheme
                input.v[i][grid_points_y-1] = -(-3*input.s[i][grid_points_y-1] + 4*input.s[i+1][grid_points_y-1] - input.s[i+2][grid_points_y-1])/(2*delta_x);
            }
            else if(i == grid_points_x-1)
            {
                //This will use the second order backward difference scheme
                input.v[i][grid_points_y-1] = -(3*input.s[i][grid_points_y-1] - 4*input.s[i-1][grid_points_y-1] + input.s[i-2][grid_points_y-1])/(2*delta_x);
            }
            else 
            {
                input.v[i][grid_points_y-1] = -(input.s[i+1][grid_points_y-1] - input.s[i-1][grid_points_y-1])/(2*delta_x);
            }
            input.u[i][grid_points_y-1] = (3*input.u[i][grid_points_y-1] - 4*input.u[i][grid_points_y-2] + input.u[i][grid_points_y-3])/(2*delta_y);
        }
        
        //Left Side
        for(int j=grid_points_y;j>=0;j--)
        {
            //X velocity will be calculate with forward difference and Y velocity will be calculated using central diffrence
            if(j == 0)
            {
                //This will use the second order forward difference scheme
                input.u[0][j] = (-3*input.s[0][j] + 4*input.s[0][j+1] - input.s[0][j+2])/(2*delta_y);
            }
            else if(j == grid_points_y-1)
            {
                //This will use the second order backward difference scheme
                input.u[0][j] = (3*input.s[0][j] - 4*input.s[0][j-1] + input.s[0][j-2])/(2*delta_y);
            }
            else 
            {
                input.u[0][j] = (input.s[0][j+1] - input.s[0][j-1])/(2*delta_y);
            }
            input.v[0][j] = - (-3*input.u[0][j] + 4*input.u[0][j] - input.u[0][j])/(2*delta_x);
        }
    */
   //Boundary values dosent have to computed since the boundary velocity are fixed.
    //Interior Points
    for(int i=1;i<grid_points_x-1;i++)
    {
        for(int j=1;j<grid_points_y-1;j++)
        {
            input.u[i][j] = (input.s[i][j+1] - input.s[i][j-1])/(2*delta_y);
            input.v[i][j] = (input.s[i+1][j] - input.s[i-1][j])/(2*delta_x);
        }
    }
    
    /*Confusions
    Now you might think that the no slip boundary coundition might run into problems by making these types of approximations with this
    definiton of velocity. but it wont happen, the reason being that if there is a no slip condition at the boundary at thge wall, then streamfunctions
    will be constant at those points, and this will result in the velcities that are again computed to be zero again at those points.*/
    return input;
}

//To compute the boundary vorticity at t2
input_stream_vorticity_2 stream_function_at_t1_to_boundary_vorticity(input_stream_vorticity_2 input,double delta_x,double delta_y)
{
    if(debugger){printf("\033[31mNow calculating boundary vorticity\033[0m\n");}
    if(debugger)
    {
        printf("\033[31mBefore updating the vorticity matrix, the stream function is as follows\033[0m\n");
        printf("\033[32m");
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0;i<grid_points_x;i++)
            {
                printf("%lf ",input.s[i][j]);
            }
            printf("\n");
        }
        printf("\n\033[0m");
    }
    //Bottom Boundary
    for(int i=0;i<grid_points_x;i++)
    {
        if(i==0) input.w[i][0] = (input.s[i][0] - 2*input.s[i+1][0] + input.s[i+2][0])/pow(delta_x,2) ;
        else if(i == grid_points_x-1) input.w[i][0] = (input.s[i-2][0] - 2*input.s[i-1][0] + input.s[i][0])/pow(delta_x,2);
        else input.w[i][0] = (input.s[i+1][0] - input.s[i-1][0])/(delta_x*2);
        input.w[i][0] = input.w[i][0] + (input.s[i][0] - 2*input.s[i][1] + input.s[i][2])/pow(delta_y,2);
        input.w[i][0] = -input.w[i][0];
    }
    if(debugger){printf("\033[31mUpdated Bottom Boundary\033[0m\n");}

    //Right Boundary
    for(int j = 0;j<grid_points_y;j++)
    {
        if(j==0) input.w[grid_points_x-1][j] = (input.s[grid_points_x-1][j] - 2*input.s[grid_points_x-1][j+1] + input.s[grid_points_x-1][j+2])/(pow(delta_y,2));
        else if(j==grid_points_y-1) input.w[grid_points_x-1][j] = (input.s[grid_points_x-1][j] - 2*input.s[grid_points_x-1][j-1] + input.s[grid_points_x-1][j-2])/pow(delta_y,2);
        else input.w[grid_points_x-1][j] = (input.s[grid_points_x-1][j+1] - 2*input.s[grid_points_x-1][j] + input.s[grid_points_x-1][j-1])/pow(delta_y,2);
        input.w[grid_points_x-1][j] = input.w[grid_points_x-1][j] + (input.s[grid_points_x-1][j] - 2*input.s[grid_points_x-2][j] + input.s[grid_points_x-3][j])/pow(delta_x,2);
        input.w[grid_points_x-1][j] = -input.w[grid_points_x-1][j];
    }
    if(debugger){printf("\033[31mUpdated Right Boundary\033[0m\n");}

    //Top Boundary (j = grid_points_y-1)
    for(int i=grid_points_x-1;i>=0;i--)
    {
        if(i==0) input.w[i][grid_points_y-1] = (input.s[i][grid_points_y-1] - 2*input.s[i+1][grid_points_y-1] + input.s[i+2][grid_points_y-1])/pow(delta_x,2);
        else if(i==grid_points_x-1) input.w[i][grid_points_y-1] = (input.s[i-2][grid_points_y-1] + 2*input.s[i-1][grid_points_y-1] + input.s[i][grid_points_y-1])/(pow(delta_x,2));
        else input.w[i][grid_points_y-1] = (input.s[i-1][grid_points_y-1] - 2*input.s[i][grid_points_y-1] + input.s[i+1][grid_points_y-1])/pow(delta_x,2);
        input.w[i][grid_points_y-1] = input.w[i][grid_points_y-1] + (input.s[i][grid_points_y-1] - 2*input.s[i][grid_points_y-2] + input.s[i][grid_points_y-3])/pow(delta_y,2);
        input.w[i][grid_points_y-1] = -input.w[i][grid_points_y-1];
    }
    if(debugger){printf("\033[31mUpdated Top Boundary\033[0m\n");}

    //left boundary
    for(int j=grid_points_y-1;j>=0;j--)
    {
        if(j==grid_points_y-1) input.w[0][j] = (input.s[0][j] - 2*input.s[0][j-1] + input.s[0][j-2])/(pow(delta_y,2));
        else if(j==0) input.w[0][j] = (input.s[0][j] - 2*input.s[0][j+1] + input.s[0][j+2])/(pow(delta_y,2));
        else input.w[0][j] = (input.s[0][j-1] - 2*input.s[0][j] + input.s[0][j+1])/(pow(delta_y,2));
        input.w[0][j] = input.w[0][j] + (input.s[0][j] - 2*input.s[1][j] + input.s[2][j])/(pow(delta_x,2));
        input.w[0][j] = -input.w[0][j];
    }
    if(debugger){printf("\033[31mUpdated Left Boundary\033[0m\n");}

    //This will print the vorticity matrix
    if(debugger)
    {
        printf("\033[31mComplete Vorticity Matrix (including boundary and interior points)\033[0m\n");
        printf("\033[32m");
        for(int j = grid_points_y-1;j >= 0;j--)
        {
            for(int i = 0;i<grid_points_x;i++)
            {
                printf("%lf ",input.w[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
    }
    return input;
}

// Step 5: The will take in input and generate the boundary condition of w at time n+1, using streamfunction at n, and solve w in the interior using the vorticity transport equation
//Scheme 1
input_stream_vorticity_2 Non_Concervative_Explicit_Upwinding_Scheme(input_stream_vorticity_2 input,double delta_x, double delta_y,double delta_t)
{
    if(debugger) printf("Into Non_Concervative_Explicit_Upwinding_Scheme\n");
    //In the function v is the dynamic viscocity
    double first_term,second_term,third_term;
    int ex,ey;
    //This will calculate the vorticity at the next time step
    for(int i = 1;i<grid_points_x-1;i++)
    {
        for(int j = 1;j<grid_points_y-1;j++)
        {
            //For each point the scheme will assume values of ex and ey. And this will depend on the velocity component at that point
            if(input.u[i][j]>0){ex=1;} else if(input.u[i][j]<0){ex=-1;} else ex=0;
            if(input.v[i][j]>0){ey=1;} else if(input.v[i][j]<0){ey=-1;} else ey=0;
            //Now we could calculate march forward in time for this point
            first_term = 0.5*(1-ex)*((input.u[i+1][j]*input.w[i+1][j] - input.u[i][j]*input.w[i][j])/delta_x) + 0.5*(1+ex)*((input.u[i][j]*input.w[i][j] - input.u[i-1][j]*input.w[i-1][j])/delta_x);
            second_term = 0.5*(1-ey)*((input.v[i][j+1]*input.w[i][j+1] - input.v[i][j]*input.w[i][j])/delta_y) + 0.5*(1+ey)*((input.v[i][j] - input.v[i][j-1]*input.w[i][j-1])/delta_y);
            third_term = dynamic_viscocity*((input.w[i+1][j] - 2*input.w[i][j]  +input.w[i-1][j])/pow(delta_x,2) + (input.w[i][j+1] - 2*input.w[i][j] + input.w[i][j-1])/pow(delta_y,2));
            input.w[i][j] = input.w[i][j] + delta_t*(third_term - first_term - second_term);
        }
    }
    return input;
}
//Scheme 2 (fix the values of b, Do not use this)
input_stream_vorticity_2 Thomas_Algorithm_SOR_by_line_Vorticity_Transport_i_intermediate(input_stream_vorticity_2 input,int k,double w,double delta_x,double delta_y,double delta_t,double re)
{ 

    int n = grid_points_x-2;
    double** A = (double**)malloc(sizeof(double*)*(n));
    
    for(int i=0;i<n;i++){A[i] = (double*)calloc(n,sizeof(double));}
    for(int row=0;row<n;row++){
        A[row][row] = 1 + (1/re)*(delta_t/(pow(delta_x,2)));
        if(row>0) A[row][row-1] = -(0.5)*(0.5*input.u[k][row]*(delta_t/delta_x) + (1/re)*(delta_t/(pow(delta_x,2))));
        else if(row<n-1) A[row][row+1] = (0.5)*(0.5*input.u[k][row+2]*(delta_t/delta_x) - (1/re)*(delta_t/(pow(delta_x,2))));
    }
    double* b = (double*)malloc(sizeof(double)*n);
    /*Changed this*/
    for(int row=0;row<n;row++){
        b[row] = (0.5)*(0.5*(input.v[k][row+1]*(delta_t/delta_y)) + (1/re)*(delta_t/pow(delta_y,2)));
        if(row == 0) b[row] = b[row] + 0.5*(0.5*input.u[k][row+1]*(delta_t/delta_x) + (1/re)*(delta_t/pow(delta_x,2)));
        else if(row == n-1) b[row] = b[row] - 0.5*(0.5*input.u[k][row+1]*(delta_t/delta_x) - (1/re)*(delta_t/pow(delta_x,2)));
    }
    
    double temp;
    for(int row=0;row<n-1;row++)
    {
        //First we make the diagonal element equal to 1
        temp= A[row][row];
        for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
        b[row] = b[row]/temp;
        
        // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
        temp = A[row+1][row];
        for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
        b[row+1] = b[row+1] - temp*b[row];        
    }
    
    input.w[k][n] = b[n-1]/A[n-1][n-1];
    for(int row=n-1;row>=1;row--)
    {
        input.w[k][row] = b[row-1]-A[row-1][row]*input.w[k][row+1];
    }
    free(A);
    free(b);
    return input;
}
input_stream_vorticity_2 Thomas_Algorithm_SOR_by_line_Vorticity_Transport_j_intermediate(input_stream_vorticity_2 input,int k,double w,double delta_x,double delta_y,double delta_t,double re)
{
    int n = grid_points_y-2;
    double** A = (double**)malloc(sizeof(double*)*(n));
    for(int i=0;i<n;i++){A[i] = (double*)calloc(n,sizeof(double));}
    
    for(int row=0;row<n;row++)
    {
        A[row][row] = 1 + (1/re)*(delta_t/(pow(delta_x,2)));
        if(row>0) A[row][row-1] = -(0.5)*(0.5*input.u[row][k]*(delta_t/delta_x) + (1/re)*(delta_t/(pow(delta_x,2))));
        else if(row<n-1) A[row][row+1] = (0.5)*(0.5*input.u[row+2][k]*(delta_t/delta_x) - (1/re)*(delta_t/(pow(delta_x,2))));
    }

    double* b = (double*)malloc(sizeof(double)*n);
    for(int row=0;row<n;row++){
        b[row] = (0.5)*(0.5*(input.v[row+1][k]*(delta_t/delta_y)) + (1/re)*(delta_t/pow(delta_y,2)));
        if(row == 0) b[row] = b[row] + 0.5*(0.5*input.u[row+1][k]*(delta_t/delta_x) + (1/re)*(delta_t/pow(delta_x,2)));
        else if(row == n-1) b[row] = b[row] - 0.5*(0.5*input.u[row+1][k]*(delta_t/delta_x) - (1/re)*(delta_t/pow(delta_x,2)));
    }
    
    double temp;
    for(int row=0;row<n-1;row++)
    {
        //First we make the diagonal element equal to 1
        temp= A[row][row];
        for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
        b[row] = b[row]/temp; 
        // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
        temp = A[row+1][row];
        for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
        b[row+1] = b[row+1] - temp*b[row];   
    }
    
    //Backward substitution
    input.w[k][n] = b[n-1]/A[n-1][n-1];
    for(int row=n-1;row>=1;row--)
    {
        input.w[k][row] = b[row-1]-A[row-1][row]*input.w[k][row+1];
    } 
    free(A);
    free(b);
    return input;
}
input_stream_vorticity_2 SOR_by_line_Vorticity_Transport_equation_intermediate(input_stream_vorticity_2 input,double delta_x,double delta_y,double delta_t)
{
    //Reynold's Number definition
    double reynold_number = (characteristic_velocity*characteristic_lenght)/dynamic_viscocity;

    if(debugger) printf("Into SOR by line function\n");
    double relaxation_factor = 0.5;

    for(int i = 1;i<grid_points_x-1;i++)
    {
        input = Thomas_Algorithm_SOR_by_line_Vorticity_Transport_i_intermediate(input,i,relaxation_factor,delta_x,delta_y,delta_t,reynold_number);
    }
    for(int j = 1;j<grid_points_y-1;j++)
    {
        input = Thomas_Algorithm_SOR_by_line_Vorticity_Transport_j_intermediate(input,j,relaxation_factor,delta_x,delta_y,delta_t,reynold_number);  
    }  

    return input;
}
//Do not use this
input_stream_vorticity_2 Thomas_Algorithm_SOR_by_line_Vorticity_Transport_i_next(input_stream_vorticity_2 input,int k,double w,double delta_x,double delta_y,double delta_t,double re)
{ 

    int n = grid_points_x-2;
    double** A = (double**)malloc(sizeof(double*)*(n));
    
    for(int i=0;i<n;i++){A[i] = (double*)calloc(n,sizeof(double));}
    for(int row=0;row<n;row++){
        A[row][row] = 1 + (1/re)*(delta_t/(pow(delta_x,2)));
        if(row>0) A[row][row-1] = (0.5)*(0.5*input.v[k][row]*(delta_t/delta_y) - (1/re)*(delta_t/(pow(delta_y,2))));
        else if(row<n-1) A[row][row+1] = -(0.5)*(0.5*input.v[k][row+2]*(delta_t/delta_y) + (1/re)*(delta_t/(pow(delta_y,2))));
    }
    double* b = (double*)malloc(sizeof(double)*n);
    /*Changed this*/
    for(int row=0;row<n;row++){
        b[row] = (0.5)*(0.5*(input.u[k][row+1]*(delta_t/delta_x)) + (1/re)*(delta_t/pow(delta_x,2)));
        if(row == 0) b[row] = b[row] + 0.5*(0.5*input.u[k][row+1]*(delta_t/delta_x) + (1/re)*(delta_t/pow(delta_x,2)));
        else if(row == n-1) b[row] = b[row] - 0.5*(0.5*input.u[k][row+1]*(delta_t/delta_x) - (1/re)*(delta_t/pow(delta_x,2)));
    }
    
    double temp;
    for(int row=0;row<n-1;row++)
    {
        //First we make the diagonal element equal to 1
        temp= A[row][row];
        for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
        b[row] = b[row]/temp;
        
        // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
        temp = A[row+1][row];
        for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
        b[row+1] = b[row+1] - temp*b[row];        
    }
    
    input.w[k][n] = b[n-1]/A[n-1][n-1];
    for(int row=n-1;row>=1;row--)
    {
        input.w[k][row] = b[row-1]-A[row-1][row]*input.w[k][row+1];
    }
    free(A);
    free(b);
    return input;
}
input_stream_vorticity_2 Thomas_Algorithm_SOR_by_line_Vorticity_Transport_j_next(input_stream_vorticity_2 input,int k,double w,double delta_x,double delta_y,double delta_t,double re)
{
    int n = grid_points_y-2;
    double** A = (double**)malloc(sizeof(double*)*(n));
    for(int i=0;i<n;i++){A[i] = (double*)calloc(n,sizeof(double));}
    
    for(int row=0;row<n;row++)
    {
        A[row][row] = 1 + (1/re)*(delta_t/(pow(delta_x,2)));
        if(row>0) A[row][row-1] = -(0.5)*(0.5*input.u[row][k]*(delta_t/delta_x) + (1/re)*(delta_t/(pow(delta_x,2))));
        else if(row<n-1) A[row][row+1] = (0.5)*(0.5*input.u[row+2][k]*(delta_t/delta_x) - (1/re)*(delta_t/(pow(delta_x,2))));
    }

    double* b = (double*)malloc(sizeof(double)*n);
    for(int row=0;row<n;row++){
        b[row] = (0.5)*(0.5*(input.v[row+1][k]*(delta_t/delta_y)) + (1/re)*(delta_t/pow(delta_y,2)));
        if(row == 0) b[row] = b[row] + 0.5*(0.5*input.u[row+1][k]*(delta_t/delta_x) + (1/re)*(delta_t/pow(delta_x,2)));
        else if(row == n-1) b[row] = b[row] - 0.5*(0.5*input.u[row+1][k]*(delta_t/delta_x) - (1/re)*(delta_t/pow(delta_x,2)));
    }
    
    double temp;
    for(int row=0;row<n-1;row++)
    {
        //First we make the diagonal element equal to 1
        temp= A[row][row];
        for(int x=0;x<n;x++){A[row][x] = A[row][x]/temp;}
        b[row] = b[row]/temp; 
        // Then we make the perform A[x][row-1] = A[x][row-1] - A[row][row-1]*A[x][row]. And this shoud be performed n-1 times
        temp = A[row+1][row];
        for(int x=0;x<n;x++){A[row+1][x] = A[row+1][x] - temp*A[row][x];}
        b[row+1] = b[row+1] - temp*b[row];   
    }
    
    //Backward substitution
    input.w[k][n] = b[n-1]/A[n-1][n-1];
    for(int row=n-1;row>=1;row--)
    {
        input.w[k][row] = b[row-1]-A[row-1][row]*input.w[k][row+1];
    } 
    free(A);
    free(b);
    return input;
}
input_stream_vorticity_2 SOR_by_line_Vorticity_Transport_equation_next(input_stream_vorticity_2 input,double delta_x,double delta_y,double delta_t)
{
    //Reynold's Number definition
    double reynold_number = (characteristic_velocity*characteristic_lenght)/dynamic_viscocity;

    if(debugger) printf("Into SOR by line function\n");
    double relaxation_factor = 0.5;

    for(int i = 1;i<grid_points_x-1;i++)
    {
        input = Thomas_Algorithm_SOR_by_line_Vorticity_Transport_i_next(input,i,relaxation_factor,delta_x,delta_y,delta_t,reynold_number);
    }
    for(int j = 1;j<grid_points_y-1;j++)
    {
        input = Thomas_Algorithm_SOR_by_line_Vorticity_Transport_j_next(input,j,relaxation_factor,delta_x,delta_y,delta_t,reynold_number);  
    }  

    return input;
}
//(Do not use this)
input_stream_vorticity_2 Non_Conervative_Implicit_ADI(input_stream_vorticity_2 input,double delta_x,double delta_y,double delta_t)
{
    //And now we will set up the thomas algorithm for internal vorticity at time n+1/2
    input = SOR_by_line_Vorticity_Transport_equation_intermediate(input,delta_x,delta_y,delta_t);
    //And this will finally advance the march to the next time n+1 
    input = SOR_by_line_Vorticity_Transport_equation_next(input,delta_x,delta_y,delta_t);
    
    return input;
}

//Scheme 3 (This will be a very long program)
//input_stream_vorticity_2 Concervative_Implicit_ADI(input_stream_vorticity_2 input);

//Main vorticity transport equation caller function
input_stream_vorticity_2 interior_vorticity_at_t1_and_boundary_vorticity_to_interior_vorticity_at_t2(input_stream_vorticity_2 input,int scheme,double delta_x,double delta_y,double delta_t)
{
    //in the function v is dynamics viscocity
    if(debugger) printf("Vorticity Transport Equation\n");
    //Non Concervative, Explicit, Upwinding Scheme
    if(scheme == 1) { input = Non_Concervative_Explicit_Upwinding_Scheme(input,delta_x,delta_y,delta_t); }
    //Non Conervative, Implicit, ADI
    else if(scheme == 2) { input = Non_Conervative_Implicit_ADI(input,delta_x,delta_y,delta_t); }
    //else if(scheme == 3) { input = Concervative_Implicit_ADI(input);}//Concervative, Implicit, ADI
    //This prints the vortocity matrix
    if(debugger)
    {
        printf("\033[31mHere is the vorticity matrix with updated interior points: \033[0m\n");
        printf("\033[32m");
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0;i<grid_points_x;i++)
            {
                printf("%lf ",input.w[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
    }
    return input;
}

//There is no need to make a new function for stream function at boundary. We could go on using the same fucntion as we did for time t1
input_stream_vorticity_2 stream_function_at_boundary_at_t2(input_stream_vorticity_2 input)
{
    return input;
}

// Step 6: This function will update the stream function equal to constant, and calculate the streamfunction at time t+1 using w at t+1 (Stream function equation)
//Again there is no need to make seperate function for this
input_stream_vorticity_2 interior_vorticity_at_t2_and_Boundary_stream_function_at_t2_to_interior_streamfunction_at_t2(input_stream_vorticity_2 input)
{

    return input;
}

// Step 7: This function will calculate the velocity and pressure at time t+1, using stream function at time t+1
input_stream_vorticity_2 streamfunction_at_t2_to_velocity_and_pressure_at_t2(input_stream_vorticity_2 input,double delta_x,double delta_y)
{
    //At the boundary we will keep the same velocity as usual, but at the interior we will update the velocity as it was
    for(int i=1;i<grid_points_x-1;i++)
    {
        for(int j=1;j<grid_points_y-1;j++)
        {
            input.u[i][j] = (input.s[i][j+1] - input.s[i][j-1])/(2*delta_y);
            input.v[i][j] = -(input.s[i+1][j] - input.s[i-1][j])/(2*delta_x);
        }
    }
    //Will calculate pressure after a while.
    //Now we will calculate the pressure matrix using the pressure poisson's equation
    if(debugger)
    {
        printf("\033[31mFinal velocity matrix for X velocity\n\033[0m");
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0;i<grid_points_x;i++)
            {
                printf("%lf\t\t",input.u[i][j]);
            }
            printf("\n");
        }
        printf("\033[31mFinal velocity matrix for Y velocity\n\033[0m");
        for(int j=grid_points_y-1;j>=0;j--)
        {
            for(int i=0;i<grid_points_x;i++)
            {
                printf("%lf\t\t",input.v[i][j]);
            }
            printf("\n");
        }
    }
    return input;
}


// This is the main function which will be used for time marching
void fluid_solver(input_stream_vorticity_2 input) 
{ 
    printf("Total Time: ");scanf("%lf",&time_max);
    printf("Time Steps: ");scanf("%d",&time_steps);
    double delta_t = time_max/time_steps;
    double delta_x = L_x/grid_points_x;
    double delta_y = L_y/grid_points_y;

    int scheme_voriticity_transport_equation;
    int scheme_streamfunction_equation;
    int scheme_pressure_poisson;
    printf("\033[31mScheme for Vorticity Transport Equation\033[0m\n"
            "1. Non Concervative, Upwinding Scheme (Explicit)\n"
            "2. Non Conervative, ADI (Implicit) \033[31m-Not in use\033[0m\n"
            "3. Concervative, ADI (Implicit) \033[31m-Not in use\033[0m\n"
    );
    printf("Enter the scheme: "); scanf("%d",&scheme_voriticity_transport_equation);

    printf("\033[31mScheme for Stream Function Equation\033[0m\n"
        "1. Sucessive Over Relaxation Method (Explicit)\n"
        "2. SOR by line (Implicit)\n"
    );
    printf("Enter the scheme: "); scanf("%d",&scheme_streamfunction_equation);

    printf("\033[31mScheme for Pressure Poisson's Equation\033[0m\n"
        "1. Sucessive Over Relaxation Method (Explicit)\n"
        "2. SOR by line (Implicit) \033[31m-Not in use\033[0m\n"
    );
    printf("Enter the scheme: "); scanf("%d",&scheme_pressure_poisson);
    printf("Enter the Dynamic Visocity : ");scanf("%lf",&dynamic_viscocity);
    printf("Enter the chracteristic lenght: "); scanf("%lf",&characteristic_lenght);
    printf("Enetr the Characteristic Velocity: ");scanf("%lf",&characteristic_velocity);

    //Here we first non dimensionalise the velocity and lenghts
    input = non_dimensionalise_velocity(input,characteristic_velocity,&delta_x,&delta_y,&delta_t,characteristic_lenght);
    //And then this velocity will be sent into the fluid solver, and all the quantites that are obtained will be non dimensionalised.

    for(int time=0;time<time_steps;time++)
    {
        if(debugger) printf("\033[1mTime step = %d\033[0m\n",time);
        
        //Since we are going to solve elliptic equation in streamfunction, it really doesnt matter if we always start with zero stream function.
        for(int i=0;i<grid_points_x;i++){for(int j=0;j<grid_points_y;j++){input.s[i][j] = 0;}}
        if(debugger) //This will print the streanm functions
        {
            printf("Here are the values of S matrix (Stream function | t = %d)\n",time);
            printf("\033[32m");
            for(int j=grid_points_y-1;j>=0;j--)
            {
                for(int i=0;i<grid_points_x;i++)
                {
                    printf("%lf ",input.s[i][j]);
                }
                printf("\n");
            }
            printf("\033[0m");
        }
        
        //This will calculate the interior values of vorticity at time t1, using velocities at time t1
        input = velocity_at_t1_to_interior_vorticity_at_t1(input,delta_x,delta_y);
        
        //This will evaluate the stream function at boundary at time t1 using velocities at the boundary
        input =  velocity_to_stream_function_at_boundary_at_t1(input,delta_x,delta_y);
        
        //This will calculate stream function at time t1 in the interior using the boundary stream function and interior vorticity
        input =  interior_vorticity_at_t1_and_boundary_stream_at_t1_to_streamfunction_at_t1(input,delta_x,delta_y,scheme_streamfunction_equation);
        
        //This will calculate recompute the velocites with this new velocities to ensure continuity and mass balance
        //input = stream_function_at_t1_to_velocity_at_t1(input,delta_x,delta_y);

        //This will compute the stream function at t1 which is completely known to calculate the vorticity at the boundary
        input =  stream_function_at_t1_to_boundary_vorticity(input,delta_x,delta_y);

        //This will compute the the vorticity at the next time step after calculating the boundary vorticity and interior vorticity at time t1
        input =  interior_vorticity_at_t1_and_boundary_vorticity_to_interior_vorticity_at_t2(input,scheme_voriticity_transport_equation,delta_x,delta_t,delta_y);
        
        //This will assign the boundary condition on the streamfunction uisng the boundary velocity to get the boundary streamfunction at tike t2
        input =  velocity_to_stream_function_at_boundary_at_t1(input,delta_x,delta_y);
        
        //This will compute the streamfunction at time t2 using the vorticity at time t2, and boundary stream fucntion at t2
        input =  interior_vorticity_at_t1_and_boundary_stream_at_t1_to_streamfunction_at_t1(input,delta_x,delta_y,scheme_streamfunction_equation);
        
        //This will genrate the pressure and velocites using the newly obtained streamfunction
        input = streamfunction_at_t2_to_velocity_and_pressure_at_t2(input,delta_x,delta_y);

        //To print the x velocity with time
        if(debugger)
        {
            printf("\033[32m");
            printf("Time Step : %d\n",time);
            for(int j=grid_points_y-1;j>=0;j--)
            {
                for(int i=0;i<grid_points_x;i++)
                {
                    printf("%lf ",input.u[i][j]);
                }
                printf("\n");
            }
            printf("\n");
            for(int j=grid_points_y-1;j>=0;j--)
            {
                for(int i=0;i<grid_points_x;i++)
                {
                    printf("%lf ",input.v[i][j]);
                }
                printf("\n");
            }
            printf("\033[0m\n");
        }
        //This is to save the file
        if(file_output)
            {
                for(int j=0;j<grid_points_y;j++)
                {
                    for(int i=0;i<grid_points_x;i++)
                    {
                        fprintf(stream_vorticity_velocity_fr,"%d,%d,%d,%lf,%lf\n",time,i,j,input.u[i][j],input.v[i][j]);
                    }
                }
            }
    }
    return;
}

/*
|**************************************************************************|
|                        PRIMITIVE VARIABLE FORMUALTION                    |
|**************************************************************************|
*/
/*
|************************** IMPORTANT COMMENTS ***********************| 
| It is important to note that only dirichlet boundary condition can  |
| be assigned on only one component of velocity.                      |
|*********************************************************************|
*/


//Important files
FILE* primitive_varibale_velocity;  //This will store all the csv data related to velocity (Complete animation)
FILE* Error_primitive_velocity;     //This will store the error between each time step 
FILE* Dilation_Error_primitive;     //Dilation error plot for each time step


//So the input will be taken into a grid, that is of the length (2*X+2)*(2*Y+2)
typedef struct input{
    double** u; //This will contain the X velocity
    double** v; //This will contain the Y velocity
    double** p; //This will contain the pressure values
    int** boundary;
    double** u_in; //This will take in the U velocity input
    double** v_in;//This will take in the V velocity input
    double** p_in; //This will be used to store the values of pressure in the final grid;
}staggered_grid;


//This function will be used to calculate the real velocity and pressure from the staggered values
void output_print(staggered_grid input){
    printf("\033[31mU_in matrix: \n");
    for(int j=grid_points_y-1;j>=0;j--)
    {
        for(int i=0;i<grid_points_x;i++)
        {
            printf("%lf ",(input.u[i][j+1] + input.u[i][j])/2);
        }
        printf("\n");
    }
    printf("\033[32mV matrix: \n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++)
        {
            printf("%lf ",(input.v[i][j] + input.v[i+1][j])/2);
        }
        printf("\n");
    }
    printf("\033[34mP matrix : \n");
    for(int j=grid_points_y-1;j>=0;j--){
        for(int i=0;i<grid_points_x;i++)
        {
            printf("%lf ",(input.p[i][j] + input.p[i+1][j] + input.p[i][j+1] + input.p[i+1][j+1])/4);
        }
        printf("\n");
    }
    printf("\033[0m");
    return;
}

//This function will be sued to non dimensialize the velocity and pressure components.
staggered_grid non_dimensionalise(staggered_grid input,double l,double v)
{
    //Non Dimensionalise U
    for(int i=0;i<grid_points_x;i++){
        for(int j=0;j<grid_points_y+1;j++){
            input.u[i][j] = input.u[i][j]/v;
        }
    }
    //Non Dimensioanlsie V
    for(int i=0;i<grid_points_x+1;i++){
        for(int j=0;j<grid_points_y;j++){
            input.v[i][j] = input.v[i][j]/v;
        }
    }
    return input;
}

//This function will give me the value of error in the divergence of velocity at each point
double dilation_error(staggered_grid input,double delta_x,double delta_y)
{
    double dilation_error = 0;
    double dilation;

    for(int i=1;i<grid_points_x-1;i++)
    {
        for(int j=1;j<grid_points_y-1;j++)
        {
            dilation = (input.u[i+1][j] - input.u[i-1][j])/(2*delta_x) + (input.v[i][j+1] - input.v[i][j-1])/(delta_y*2);
            dilation_error = dilation_error + dilation;
        }
    }
    return dilation_error;
}

//This function will solve the momentum eqaution
staggered_grid momentum(staggered_grid input,double delta_x,double delta_y,double delta_t,double reynolds)
{
    if(debugger) printf("Into the Momentum Function\n");
    double dell_velocity_squared;
    double dell_velocity_product;
    double dell_pressure;
    double double_dell_velocity_x,double_dell_velocity_y;
    //So this will solve the velocity at the interior points only (X velocity)
    //And since the velocity in x is define for the matrix if dimension (Nx)*(Ny+1) therefore X = [1,Nx-1], and Y = [1,Ny]
    //We will first take the duplicate of U
    double** u_old = (double**)malloc(sizeof(double*)*grid_points_x);
    for(int i=0;i<grid_points_x;i++){u_old[i] = (double*)malloc(sizeof(double)*(grid_points_y+1));}
    for(int i=0;i<grid_points_x;i++){for(int j=0;j<grid_points_y+1;j++){u_old[i][j] = input.u[i][j];}}
    
    if(strong_debugger) printf("Computing Interior U matrix\n");
    //Solver (X velocity)
    for(int i=1;i<grid_points_x-1;i++)
    {
        for(int j=1;j<grid_points_y;j++)
        {
            if(strong_debugger) printf("i = %d | j = %d\n",i,j);
            dell_velocity_squared = ((pow(input.u[i+1][j] + input.u[i][j],2) - pow(input.u[i][j] + input.u[i-1][j],2))/4*delta_x);
            dell_velocity_product = (((input.u[i][j] + input.u[i][j+1])/2)*((input.v[i][j] + input.v[i+1][j])/2) - ((input.u[i][j] + input.u[i][j-1])/2)*((input.v[i][j-1] + input.v[i+1][j-1])/2))/delta_y;                                 
            dell_pressure = (input.p[i+1][j] - input.p[i][j])/delta_x;
            double_dell_velocity_x = (input.u[i+1][j] + input.u[i-1][j] - 2*input.u[i][j])/(pow(delta_x,2));
            double_dell_velocity_y = (input.u[i][j+1] + input.u[i][j-1] - 2*input.u[i][j])/(pow(delta_y,2));
        
            //now we update the x velocity at the interior
            input.u[i][j] = input.u[i][j] + delta_t*(-dell_pressure - dell_velocity_squared - dell_velocity_product + (1/reynolds)*(double_dell_velocity_x + double_dell_velocity_y));
            if(strong_debugger) printf("u[%d][%d] = %lf\n",i,j,input.u[i][j]);
        }
    }
    //And now we will set the boundary values of U
    for(int i=0;i<grid_points_x;i++)
    { 
        //Bottom boundary
        input.u[i][0] = 2*input.u_in[i][0]-input.u[i][1];
        //Top boundary
        input.u[i][grid_points_y] = 2*input.u_in[i][grid_points_y-1] -input.u[i][grid_points_y-1];
    }
    for(int j=0;j<grid_points_y;j++)
    {
        //Left Boundary
        if(j==grid_points_y) input.u[0][j] = 2*input.u_in[0][j];
        else input.u[0][j] = 2*input.u_in[0][j] - input.u[0][j+1];
        //Right
        if(j==grid_points_y) input.u[grid_points_x-1][j] = 2*input.u_in[grid_points_x-1][j];
        else input.u[grid_points_x-1][j] = 2*input.u_in[grid_points_x-1][j] - input.u[grid_points_x-1][j+1];
    }

    //And now we will update the y velocity at the interior
    if(strong_debugger) printf("Computing Interior V matrix\n");
    if(strong_debugger)
    {
        printf("We will first print the entire V velocity matrix\n");
        for(int i=0;i<grid_points_x+1;i++){for(int j=0;j<grid_points_y;j++){printf("v[%d][%d] = %lf\n",i,j,input.v[i][j]);}}
    }
    //Solver (Y velocity)
    for(int i=1;i<grid_points_x;i++){
        for(int j=1;j<grid_points_y-1;j++){
            if(strong_debugger) printf("i = %d | j = %d\n",i,j);
            dell_velocity_squared = (pow((input.v[i][j+1] + input.v[i][j]),2) - pow((input.v[i][j] + input.v[i][j-1]),2))/(4*delta_y);
            dell_velocity_product = ((u_old[i][j+1] + u_old[i][j])*(input.v[i][j] + input.v[i+1][j]) - (u_old[i-1][j+1] + u_old[i-1][j])*(input.v[i][j] + input.v[i-1][j]))/(4*delta_x);    
            dell_pressure = (input.p[i][j+1] - input.p[i][j])/delta_y;
            double_dell_velocity_x = (input.v[i+1][j] + input.v[i-1][j] - 2*input.v[i][j])/(pow(delta_x,2));
            double_dell_velocity_y = (input.v[i][j+1] + input.v[i][j-1] - 2*input.v[i][j])/(pow(delta_y,2));

            //now we update the y velocity at the interior
            input.v[i][j] = input.v[i][j] + delta_t*(-dell_velocity_product - dell_velocity_squared - dell_pressure + (1/reynolds)*(double_dell_velocity_x + double_dell_velocity_y));
            if(strong_debugger) printf("v[%d][%d] = %lf\n",i,j,input.v[i][j]);
        }
    }
    
    //And now we will set the boundary values on V
    for(int i=0;i<grid_points_x;i++)
    {
        //Top boundary
        if(i==0) input.v[i][grid_points_y-1] = 2*input.v_in[i][grid_points_y-1];
        else input.v[i+1][grid_points_y-1] = 2*input.v_in[i][grid_points_y-1] - input.v[i][grid_points_y-1];
        //Bottom Boundary
        if(i==0) input.v[i][0] = 2*input.v_in[i][0];
        else input.v[i+1][0] = 2*input.v_in[i][0] - input.v[i][0];
    }
    for(int j=0;j<grid_points_y;j++)
    {
        //Right Boundary
        input.v[0][j] = 2*input.v_in[0][j] -input.v[1][j];
        //Left Boundary
        input.v[grid_points_x][j] = 2*input.v_in[grid_points_x-1][j] -input.v[grid_points_x-1][j];
    }
    
    if(debugger) printf("Out of Momentum Function\n");
    //Free the u_old array
    for(int i=0;i<grid_points_x;i++){free(u_old[i]);} free(u_old);
    return input;
}
//This function will take in the input and will solve pressure at time n+1
staggered_grid pressure(staggered_grid input,double delta_x,double delta_y,double delta_t,double rho)
{
    if(debugger) printf("Into pressure function\n");
    double dell_u;
    double dell_v;
    for(int i=1;i<grid_points_x;i++)
    {
        for(int j=1;j<grid_points_y;j++)
        {
            dell_u = (input.u[i][j] - input.u[i-1][j])/delta_x;
            dell_v = (input.v[i][j] - input.v[i][j-1])/delta_y;
            input.p[i][j] = input.p[i][j] - rho*delta_t*(dell_u + dell_v);
        }
    }
    //And now we will set the boundary condition
    for(int i=0;i<grid_points_x+1;i++)
    {
        input.p[i][0] = input.p[i][1];
        input.p[i][grid_points_y] = input.p[i][grid_points_y-1];
    }
    for(int j=0;j<grid_points_y+1;j++)
    {
        input.p[0][j] = input.p[1][j];
        input.p[grid_points_x][j] = input.p[grid_points_x-1][j];
    }
    if(debugger) printf("Out of pressure function\n");
    return input;
}

//This function will make U,V from U_in and V_in
staggered_grid stagger_maker(staggered_grid input)
{
    /*
    FILE* u_in_matrix;
    FILE* v_in_matrix;
    u_in_matrix = fopen("U_in_matrix.txt","w");
    v_in_matrix = fopen("V_in_matrix.txt","w");
    for(int j=grid_points_y-1;j>=0;j--)
    {
        for(int i=0;i<grid_points_x;i++)
        {
            fprintf(u_in_matrix,"%lf ",input.u_in[i][j]);
            fprintf(v_in_matrix,"%lf ",input.v_in[i][j]);
        }
        fprintf(u_in_matrix,"\n");
        fprintf(v_in_matrix,"\n");
    }
    fclose(u_in_matrix);
    fclose(v_in_matrix);
    */

    //Now we will apply the algorithm for filling the values of u and v from u_in and v_in
    for(int i=0;i<grid_points_x;i++)
    {
        // We first update the top points
        input.u[i][grid_points_y] = 2*input.u_in[i][grid_points_y-1];
        // And then we propogate downwards
        for(int j=grid_points_y-1;j>=0;j--){input.u[i][j] = 2*input.u_in[i][j] - input.u[i][j+1];}
    }
    for(int j=0;j<grid_points_y;j++)
    {
        input.v[0][j] = 2*input.v_in[0][j];
        for(int i=0;i<grid_points_x;i++){input.v[i+1][j] = 2*input.v_in[i][j] - input.v[i][j];}
    }
    
    /*
    //We may print the values of velocity and pressure.
    int ask;
    printf("Do you want to print the values of velocity and pressure (0 No | 1 Yes)?");
    scanf("%d",&ask);
    if(ask)
    {
        //Now we will print the velocity and pressure matrix 
        printf("The Staggered matrix are named as U and V and P\n"
            "And Non-Staggered gird is marked as U_in and V_in and P_in\n");
        printf("\033[31mU Matrix : \n");
        for(int j = grid_points_y;j>=0;j--)
        {
            for(int i=0;i<grid_points_x;i++)
            {
                printf("%lf ",input.u[i][j]);
            }
            printf("\n");
        }
        printf("\033[32mV Matrix : \n");
        for(int j = grid_points_y-1;j>=0;j--){
            for(int i=0;i<grid_points_x+1;i++)
            {
                printf("%lf ",input.v[i][j]);
            }
            printf("\n");
        }
        printf("\033[34mP Matrix : \n");
        for(int j = grid_points_y;j>=0;j--){
            for(int i=0;i<grid_points_x+1;i++)
            {
                printf("%lf ",input.p[i][j]);
            }
            printf("\n");
        }
        printf("\033[0m");
    }
    */
    return input;
}

//This is the main solver for the staggered grid method, and it will be take in the output_staggered_grid type datatype and will output the values in the form of csv data
void staggered_grid_fluid_solver(staggered_grid input,bool isTest,double U1,double V1,double U2,double V2,double U3,double V3)
{
    printf("Dimension of domian in X : "); fflush(stdout); scanf("%lf",&L_x);
    printf("Dimension of domain in Y : "); fflush(stdout); scanf("%lf",&L_y);
    printf("Total Time : "); fflush(stdout); scanf("%lf",&time_max);

    double delta_x = L_x/(grid_points_x-1);
    double delta_y = L_y/(grid_points_y-1);
    double delta_t = time_max/time_steps;

    double characteristic_lenght; printf("Characteristic Length : "); fflush(stdout); scanf("%lf",&characteristic_lenght);
    double characteristic_velocity; printf("Characteristic Velocity : "); fflush(stdout); scanf("%lf",&characteristic_velocity);
    double dynamic_viscocity; printf("Dynamic viscocity : "); fflush(stdout); scanf("%lf",&dynamic_viscocity);
    double reynold = (characteristic_velocity*characteristic_lenght)/dynamic_viscocity;

    input = non_dimensionalise(input,characteristic_lenght,characteristic_velocity);
    double rho; printf("Whats the artificial compressiblity : ");fflush(stdout); scanf(" %lf",&rho);

    //Initalise the final pressure matrix
    input.p_in  = (double**)malloc(sizeof(double*)*grid_points_x);
    for(int i=0;i<grid_points_x;i++){input.p_in[i] = (double*)calloc(sizeof(double),grid_points_y);}
    
    //Used for calculating the error
    double** velocity_old = (double**)malloc(sizeof(double)*(grid_points_x));
    for(int i=0;i<grid_points_x;i++){velocity_old[i] = (double*)malloc(sizeof(double)*grid_points_y);}
    double* momentum_error = (double*)malloc(sizeof(double)*time_steps); //This will store the error at each point
    double* continuity_error = (double*)malloc(sizeof(double)*time_steps); //this will be sued to keep the track of continuity erro at each time step

    //For error temp variables
    double u_velocity;
    double v_velocity;

    //This will be used to check if the ouptut has blown up
    bool blown = false;

    time_t start = clock(); //This will start the clock

    //Time Marching
    for(int time=0;time<time_steps;time++)
    {
        //We will first store all the values of velocity in velocity old
        for(int i=0;i<grid_points_x;i++)
        {
            for(int j=0;j<grid_points_y;j++)
            {
                u_velocity = (input.u[i][j] + input.u[i][j+1])/2;
                v_velocity = (input.v[i][j] + input.v[i+1][j])/2;
                velocity_old[i][j] = pow((pow(u_velocity,2) + pow(v_velocity,2)),0.5);
            }
        }
        
        if(time%1==0){ printf("\033[32mTime = %d\n\033[0m",time);}
        
        
        //Then we recompute the u,v from u_in and v_in

        //First we solve the momentum eqaution of the interior velocity
        input = momentum(input,delta_x,delta_y,delta_t,reynold);
        
        //Then we solve the continuity eqaution for pressure
        input = pressure(input,delta_x,delta_y,delta_t,rho);

        printf("Instability Checker\n");
        //Before proceeding to the next time step we will check if all the enteries inside the velocity is valid number and not eqaul to nan
        for(int i=0;i<grid_points_x;i++){
            for(int j=0;j<grid_points_y;j++){
                double x_vel = (input.u[i][j] + input.u[i][j+1])/2;
                double y_vel = (input.v[i][j] + input.v[i+1][j])/2;
                if(isnan(x_vel) || isnan(y_vel)){
                    blown = true;
                    break;
                }
            }
        }
        if(blown) break;

        //And then we calculate the error after iterations (For steady state)
        for(int i=0;i<grid_points_x;i++)
        {
            for(int j=0;j<grid_points_y;j++)
            {
                u_velocity = (input.u[i][j] + input.u[i][j+1])/2;
                v_velocity = (input.v[i][j] + input.v[i+1][j])/2;
                momentum_error[time] = momentum_error[time] + pow((velocity_old[i][j] - pow((pow(u_velocity,2) + pow(v_velocity,2)),0.5)),2);
            }
        }
        printf("\033[31mMomentum_Error = %lf\n\033[0m",momentum_error[time]);
        //calculate the continuity error at each time
        continuity_error[time] = dilation_error(input,delta_x,delta_y);
        printf("\033[31mContinuity_Error = %lf\n\033[0m",continuity_error[time]);

        //This will store the value of velocity, dilation and momentum error in file 'primitive_varibale_velocity','Error_primitive_velocity','Conti_Error_primitive_velocity'
        if(file_output){
            for(int i=0;i<grid_points_x;i++){
                for(int j=0;j<grid_points_y;j++){
                    fprintf(primitive_varibale_velocity,"%d,%d,%d,%lf,%lf\n",time,i,j,(input.u[i][j+1] + input.u[i][j])/2,(input.v[i][j] + input.v[i+1][j])/2);
                    fprintf(Error_primitive_velocity,"%d,%lf\n",time,momentum_error[time]);
                    fprintf(Dilation_Error_primitive,"%d,%lf\n",time,continuity_error[time]);
                }
            }
        }

        //We could also update the values of u_in and v_in to the new ones
        for(int i=0;i<grid_points_x;i++){
            for(int j=0;j<grid_points_y;j++){
                input.u_in[i][j] = (input.u[i][j+1] + input.u[i][j])/2;
                input.v_in[i][j] = (input.v[i][j] + input.v[i+1][j])/2;
                input.p_in[i][j] = (input.p[i][j] + input.p[i+1][j] + input.p[i][j+1] + input.p[i+1][j+1])/4;
            }
        }

        //Then we re-enfore the boundary condition
        for(int i=0;i<grid_points_x;i++)
        {
            for(int j=0;j<grid_points_y;j++)
            {
                if(input.boundary[i][j] == 1){
                    input.u_in[i][j] = U1;
                    input.v_in[i][j] = V1;
                }
                else if(input.boundary[i][j] == 2)
                {
                    input.u_in[i][j] = U2;
                    input.v_in[i][j] = V2;
                }
                else if(input.boundary[i][j] == 3)
                {
                    input.u_in[i][j] = U3;
                    input.v_in[i][j] = V3;
                }
            }
        }

        //After enforcing the boundary condition you again make the U,V from the new U_in and v_in
        input = stagger_maker(input);

        if(debugger)
        {
            //Final values could be printed
            printf("\033[31mX-Matrix\n");
            for(int j=grid_points_y-1;j>=0;j--){for(int i=0;i<grid_points_x;i++){printf("%lf ",input.u_in[i][j]);} printf("\n");}
            printf("\033[32mY-Matrix\n");
            for(int j=grid_points_y-1;j>=0;j--){for(int i=0;i<grid_points_x;i++){printf("%lf ",input.v_in[i][j]);} printf("\n");}
            printf("\033[34mP-Matrix\n");
            for(int j=grid_points_y-1;j>=0;j--){for(int i=0;i<grid_points_x;i++){printf("%lf ",input.p_in[i][j]);} printf("\n");}
        }

        //Next iteration...
    }

    if(!blown) //if the solution is not blown up
    {
        time_t end = clock();
        if(computation_time){
            printf("\033[34m");
            double cpu_time = (double)(end-start)/CLOCKS_PER_SEC;
            printf("CPU Time = %lf\n",cpu_time);
            printf("\033[0m");
        }
        if(plot){system("python Error_Primitve.py");}
    }
    else if(blown) printf("Solution has blown up\n"); //If the solution is blown up

    for(int i=0;i<grid_points_x;i++){free(velocity_old[i]);} free(velocity_old); //This is to free the old velocity
    free(momentum_error); //This is to free the error.
    free(continuity_error);

    return;
}

//Post processing Function
void post_processing(staggered_grid input)
{

}

//This function will take the input from the user about the problem
void primitive_variable_formulation()
{
    printf("Primitve Variable Formulation function\n");
    //Grid points size
    printf("Size of grid in X : "); fflush(stdout); scanf("%d",&grid_points_x);
    printf("Size of grid in Y : "); fflush(stdout); scanf("%d",&grid_points_y);
    printf("Number of time steps : "); fflush(stdout); scanf("%d",&time_steps);

    //Input is initalised here
    staggered_grid input;
    //In Staggered grid method, we have the row is u matirx to be (Nx+1) and columns to be (Ny)
    input.u = (double**)malloc(sizeof(double)*(grid_points_x));for(int i =0;i<(grid_points_x);i++){input.u[i] = (double*)malloc(sizeof(double)*(grid_points_y+1));}
    //In staggered grid method, we have the columns in V matrix to be (Ny+1) amd row to be (Nx)
    input.v = (double**)malloc(sizeof(double)*(grid_points_x+1));for(int i =0;i<grid_points_x+1;i++){input.v[i] = (double*)malloc(sizeof(double)*(grid_points_y));}
    //And pressure will be of the size (Nx+1) times (Ny+1)
    input.p = (double**)malloc(sizeof(double)*(grid_points_x+1));for(int i =0;i<grid_points_x+1;i++){input.p[i] = (double*)malloc(sizeof(double)*(grid_points_y+1));}

    //Creates the input matrix
    input.u_in = (double**)malloc(sizeof(double*)*grid_points_x);
    input.v_in = (double**)malloc(sizeof(double*)*grid_points_x);
    for(int i=0;i<grid_points_x;i++)
    {
        input.u_in[i] = (double*)malloc(sizeof(double)*grid_points_y);
        input.v_in[i] = (double*)malloc(sizeof(double)*grid_points_y);
    }
    
    //And the pressure will always be initialised by 1
    for(int i=0;i<grid_points_x+1;i++){for(int j=0;j<grid_points_y+1;j++){input.p[i][j] = 1;}}
    
    //And now we could start filling the values from the user. But here we give the option of whether the user want to fill the values manually or wnats to use the test problems
    const char* options[] = {"1.Test Problem","2.Custom Problem (Manual)","3.Custom Problem (Image)"};
    int choice = 0;
    int total_choices = 3;
    while(1){
        system("cls");
        printf("\033[1;31m\tVersion Selection\033[0m\n");
        printf("Use arrow keys and Enter to select:\n");
        for(int i=0;i<total_choices;i++)
        {
            if (i == choice) printf("\t\t\033[1;32m> %s <\033[0m\n", options[i]);
            else printf("\t\t\033[31m %s\033[0m\n", options[i]);
        }

        int c_version = _getch();
        if (c_version == 0 || c_version == 224) {  // arrow keys are 2-char codes
            switch (_getch()) {
                case 72: choice = (choice - 1 + total_choices) % total_choices; break; // up
                case 80: choice = (choice + 1) % total_choices; break;     // down
            }
        } else if (c_version == 13) {  // Enter key
            printf("\nYou selected:\033[1m %s\n\033[0m", options[choice]);
            break;
        }
    }
    //The choosen choice is
    switch (choice+1){
        /*Custom Problem (Image)*/
        case 3:{
            //This will call the python command, and the code will prepare the boundary matrix for later use.
            char command[100];
            sprintf(command,"python Image_Processing/input.py %d %d",grid_points_x,grid_points_y);
            system(command);

            //Memeory Allocation to the boundary matrix.
            input.boundary = (int**)malloc(sizeof(int*)*grid_points_x);
            for(int i=0;i<grid_points_x;i++){input.boundary[i] = (int*)malloc(sizeof(int)*grid_points_y);}
            //Taking values from the boundary matrix.
            FILE* boundary_file = fopen("./Image_Processing/boundary_matrix.txt","r");
            if(boundary_file == NULL){printf("Error opening boundary input matrix\n");}
            else{
                printf("Now the input is being scanned\n");
                for(int j=grid_points_y-1;j>=0;j--)
                {
                    for(int i=0;i<grid_points_x;i++)
                    {
                        fscanf(boundary_file,"%d",&input.boundary[i][j]);
                    }
                }
            }

            //Now we ask for the boundary values of the velocity
            printf("Important to note that the White is assumed to be interior point\n");
            double U_1,V_1;//Dirichlet boundary condition on the flow with boundary value 1
            double U_2,V_2;//Dirichlet boundary condition on the flow with boundary value 2
            double U_3,V_3;//Dirichlet boundary condition on the flow with boundary value 3
            printf("Dirichlet condition on the Black color (U V): ");scanf("%lf %lf",&U_1,&V_1);
            printf("Dirichlet condition on the Blue color (U V): ");scanf("%lf %lf",&U_2,&V_2);
            printf("Dirichlet condition on the Red color (U V): ");scanf("%lf %lf",&U_3,&V_3);

            //Preparation of U_in and V_in matrix
            for(int i=0;i<grid_points_x;i++)
            {
                for(int j=0;j<grid_points_y;j++)
                {
                    //printf("boundary[%d][%d] = %d\n",i,j,input.boundary[i][j]);
                    //If the value in the boundary matrix is equal to 0
                    if(input.boundary[i][j] == 0)
                    {
                        input.u_in[i][j] = 0;
                        input.v_in[i][j] = 0;
                    }
                    else if(input.boundary[i][j] == 1)
                    {
                        input.u_in[i][j] = U_1;
                        input.v_in[i][j] = V_1;
                    }
                    else if(input.boundary[i][j] == 2)
                    {
                        input.u_in[i][j] = U_2;
                        input.v_in[i][j] = V_2;
                    }
                    else if(input.boundary[i][j] == 3)
                    {
                        input.u_in[i][j] = U_3;
                        input.v_in[i][j] = V_3;
                    }
                }
            }

            //To prepare the U, and V matrix
            printf("Into the staggered grid maker function\n");
            input = stagger_maker(input);

            //After the staggered grid is taken into the u_in, and v_in and fta
            staggered_grid_fluid_solver(input,false,U_1,V_1,U_2,V_2,U_3,V_3);
            break;
        }
        /*Custom Problem (Manual)*/
        case 2:
        {
            printf("Enter the U matrix :\n");//Taking in the U matirx
            for(int j=grid_points_y-1;j>=0;j--){
                for(int i=0;i<grid_points_x;i++){
                    scanf("%lf",&input.u_in[i][j]);
                    //if(j==grid_points_y-1) input.u_in[i][j] = 10;
                    //else input.u_in[i][j] = 0;
                }
            }
            printf("Enter the V matrix :\n");//Taking in the V matirx
            for(int j=grid_points_y-1;j>=0;j--){
                for(int i=0;i<grid_points_x;i++){
                    scanf("%lf",&input.v_in[i][j]);
                    //if(i==0) input.v_in[i][j] = 10;
                    //else input.v_in[i][j] = 0;
                }
            }
            
            //For making the staggered grid values
            input = stagger_maker(input);
            output_print(input);

            printf("Into the Fluid Solver\n");
            //And now this input could be sent into fluid fluid_solver
            staggered_grid_fluid_solver(input,false,0,0,0,0,0,0);

            /*Here we don't really need the dirchlet velocity at all*/
            break;
        }
        /*Test Problem */
        case 1:
            //Now we set the input matrix
            for(int i=0;i<grid_points_x;i++){for(int j=0;j<grid_points_y+1;j++){if(j==grid_points_y){input.u[i][j] = 20;}else input.u[i][j] = 0;}}
            for(int i=0;i<grid_points_x+1;i++){for(int j=0;j<grid_points_y;j++)input.v[i][j]=0;}
            for(int i=0;i<grid_points_x+1;i++){for(int j=0;j<grid_points_y+1;j++)input.p[i][j]=1;}
            //And now this input could be sent into fluid fluid_solver
            staggered_grid_fluid_solver(input,false,0,0,0,0,0,0);
            break;
        default:
            break;
    }

    //Free the u,v,p
    for(int i=0;i<grid_points_x;i++){free(input.u[i]);}
    for(int i=0;i<grid_points_x+1;i++){free(input.v[i]); 
    free(input.p[i]);}
    free(input.u);free(input.v);free(input.p);

    //free u_in, v_in and p_in
    for(int i=0;i<grid_points_x;i++){free(input.u_in[i]);free(input.v_in[i]);free(input.p_in[i]);}
    free(input.u_in);free(input.v_in);free(input.p_in); 
    
    //Free boundary
    for(int i=0;i<grid_points_x;i++){free(input.boundary[i]);}
    free(input.boundary);
}

/* MAIN FUNCTION ************************************************************************************/
int main(){ 
    printf("\033[0m");
    system("cls");
    opening_animation_2();
    system("cls");

    const char *version[] = {"Version 1 (Toy Problems)","Version 2 (Stream Vorticity Formulation)","Version 3 (Stagered Grid Formulation)"};
    int choice_version = 0;
    int version_number = 3;
    while(1){
        system("cls");
        printf("\033[1;31m\tVersion Selection\033[0m\n");
        printf("Use arrow keys and Enter to select:\n");
        for(int i=0;i<version_number;i++)
        {
            if (i == choice_version) printf("\t\t\033[1;32m> %s <\033[0m\n", version[i]);
            else printf("\t\t\033[31m %s\033[0m\n", version[i]);
        }

        int c_version = _getch();
        if (c_version == 0 || c_version == 224) {  // arrow keys are 2-char codes
            switch (_getch()) {
                case 72: choice_version = (choice_version - 1 + version_number) % version_number; break; // up
                case 80: choice_version = (choice_version + 1) % version_number; break;     // down
            }
        } else if (c_version == 13) {  // Enter key
            printf("\nYou selected:\033[1m %s\n\033[0m", version[choice_version]);
            break;
        }
    }

    switch (choice_version+1)
    {
        //Toy Problems
        case 1:{
            input_type_of_problem();
            if(dimension==1){
                input_1 input;
                input = input_data_1(input);
                output_1(input);
                free(input.boundary_matrix);
                free(input.value_matrix);
            }
            else if(dimension==2)
            {
                input_2 input;
                input = input_data_2(input);
                output_2(input);
                free(input.boundary_matrix);
                free(input.value_matrix);
            }
            else if(dimension==3) printf("Computation for 3D problem is not supported yet\n");
            else{
                printf("Enter valid dimension\n");
                main();
                return 0;
            }
            break;
        }
        //Stream Vorticity Formulation
        case 2:{
            if(file_output){
                stream_vorticity_velocity_fr = fopen("D:/CloseFOAM/Stream_Function_Vorticity_formulation_Results_Velocity.csv","w");
                
                if(stream_vorticity_velocity_fr == NULL)
                {
                    printf("The file has not loaded properly.\n"
                        "The output genrated my not be saved in the file\n");
                }
                else printf("File will be generated\n");
            }
            time_t program_start = clock(); //Starts the clock
            input_stream_vorticity_2 input; //Initialise the input
            
            //Option to take the input
            const char *options[] = {"Custom_Terminal_Input", "Test_Input(Cavity_lid_problem)","File_input (Do not use)"};
            int choice = 0;
            int n = 3;
            while (1) {
                system("cls"); // clear screen
                logo();
                printf("\033[1;34m\nVersion %d : %s\033[0m\n",choice_version+1,version[choice_version]);
                printf("\nUse arrow keys and Enter to select:\n");

                for (int i = 0; i < n; i++) {
                    if (i == choice) printf("\033[1;32m> %s <\033[0m\n", options[i]);
                    else printf("\033[2;31m %s\033[0m\n", options[i]);
                }

                int c = _getch();
                if (c == 0 || c == 224) {  // arrow keys are 2-char codes
                    switch (_getch()) {
                        case 72: choice = (choice - 1 + n) % n; break; // up
                        case 80: choice = (choice + 1) % n; break;     // down
                    }
                } else if (c == 13) {  // Enter key
                    printf("\nYou selected:\033[1m %s\n\033[0m", options[choice]);
                    break;
                }
            }
            switch (choice+1)
            {
                case 1:
                    input = input_data_stream_function_vorticity_2(input);
                    break;
                case 2:
                    input = test_input_function(input,version,options,choice_version,choice);
                case 3:
                    //printf("In Progress\n");
                    break; 
                default:
                    break;
            }
            //This is to proceed solving the equations
            if(!terminate) fluid_solver(input);
            if(terminate) {printf("\033[1;31mThe program was terminated\033[0m\n"); return 0;}
            
            time_t program_end = clock(); //Ends the clock
            if(computation_time) //Prints the total CPU time
            {
                printf("\n\033[34mTotal compuation time\n");
                double total_program_time = (double)(program_end-program_start)/CLOCKS_PER_SEC;
                printf("Total Time Taken : %lf",total_program_time);
                printf("\n\033[0m");
            }
            
            //This will run the python script to visulaise the animation of velocity. The frames will also be saved in the folder
            if(plot){
                char command[256];
                sprintf(command, "cd /d D:\\CloseFOAM && python Velocity_Script.py %d %d", grid_points_x, grid_points_y);
                int ret = system(command);
                if(ret != 0) {printf("Error running Python script!\n");}
            }
            //It is important to free the input
            if(&input!=NULL){
                for(int i=0;i<grid_points_x;i++)
                {
                    free(input.boundary_matrix[i]);
                    free(input.p[i]);
                    free(input.s[i]);
                    free(input.u[i]);
                    free(input.v[i]);
                    free(input.w[i]);
                }
                free(input.boundary_matrix);
                free(input.p);
                free(input.s);
                free(input.u);
                free(input.v);
                free(input.w);
            }
            if(file_output) {fclose(stream_vorticity_velocity_fr);}//Closing file
            
            break;
        }
        //Staggered Grid Formulation
        case 3:{
            //Opens the file
            if(file_output)
            {
                primitive_varibale_velocity = fopen("D:/CloseFOAM/Primitive_varibale_Results_Velocity.csv","w");
                Error_primitive_velocity = fopen("D:/CloseFOAM/Primitive_varibale_Error.csv","w");
                Dilation_Error_primitive = fopen("D:/CloseFOAM/Primitive_varibale_Dilation_Error.csv","w");
                if(primitive_varibale_velocity == NULL || Error_primitive_velocity == NULL || Dilation_Error_primitive == NULL)
                {
                    printf("The file has not loaded properly.\n"
                        "The output genrated may not be saved in the file\n");
                }
                else printf("File will be generated\n");
            }

            //Solver
            primitive_variable_formulation();

            //This will plot the animation
            if(plot){
                printf("About to plot the result\n");
                char command[256];
                sprintf(command, "python Python_Scripts/Velocity_Script_primitive.py %d %d", grid_points_x, grid_points_y);
                int ret = system(command);
                if(ret != 0) {printf("Error running Python script!\n");}
            }

            //Closing the files
            if(file_output) {
                fclose(primitive_varibale_velocity);
                fclose(Error_primitive_velocity);
            }
        }
        default:
            break;     
    }
    printf("\n\033[1;31mEnd\033[0m\n");
    return 0;
}
