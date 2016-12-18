!
! Abstract Data Type Stack Implemented in Fortran.
! (c) 2004-2010, Athanasios Migdalas @ AUTh
! version 0.5 @ 2010/05/15

MODULE StackObject
    USE global
  IMPLICIT NONE

  !-------------! Interface provided to clients !-------------------!

  TYPE StackT ! Data Type Stack, 
     PRIVATE  ! all members private and inaccessible to clients
     TYPE(solution), ALLOCATABLE, DIMENSION(:) :: Key
     INTEGER :: Last, ErrorNr, Size, Count
     LOGICAL :: Error
  END TYPE StackT

  ! Stack Operators that are available to clients under different names

  INTERFACE push
     MODULE PROCEDURE add_to_stack
  END INTERFACE

  INTERFACE pop
     MODULE PROCEDURE delete_from_stack
  END INTERFACE

  INTERFACE top
     MODULE PROCEDURE select_from_stack
  END INTERFACE
  
  !INTERFACE scount
  !   MODULE PROCEDURE getCount
  !END INTERFACE

  ! These Operators are not available to clients

  PRIVATE :: add_to_stack, delete_from_stack, select_from_stack, SetErrorFlags

  ! All other Operators are public and thus, available to clients

CONTAINS  !---------------! Operators are implemented below !---------------!

  ! (1) Construct and initialize a new Stack

  FUNCTION NewStack(StackSize) RESULT(Stack)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: StackSize
    INTEGER :: ierror
    TYPE(StackT) :: Stack

    IF ( .NOT. ALLOCATED(Stack%Key) ) THEN
       ALLOCATE( Stack%Key(StackSize), STAT=ierror)
       IF ( ierror /= 0 ) THEN
          CALL SetErrorFlags(Stack,.TRUE.,-1)
       ELSE
          !Stack%Key = 0
          Stack%Size = StackSize
          Stack%Last = 0
          Stack%Count = 0
          CALL SetErrorFlags(Stack,.FALSE.,0)
       ENDIF
    ENDIF

  END FUNCTION NewStack

  ! (2) Destroy a Stack and free the memory it occupies

  SUBROUTINE DisposeStack(Stack)
    IMPLICIT NONE
    INTEGER :: ierror
    TYPE(StackT) :: Stack

    CALL SetErrorFlags(Stack,.FALSE.,0)

    IF ( ALLOCATED(Stack%Key) ) THEN
       DEALLOCATE( Stack%Key, STAT=ierror)
       IF ( ierror /= 0 ) THEN
          CALL SetErrorFlags(Stack,.TRUE.,-2)
       ENDIF
    ENDIF

  END SUBROUTINE DisposeStack

  ! (3) Add a Node to the Stack (Push)

  SUBROUTINE add_to_stack(Node,Stack)
    IMPLICIT NONE
    TYPE(solution) :: Node
    TYPE(StackT) :: Stack

    IF ( Stack%Last >= Stack%Size ) THEN
       CALL SetErrorFlags(Stack,.TRUE.,-3)
    ELSE
       Stack%Last = Stack%Last + 1
       Stack%Key( Stack%Last ) = Node
       Stack%Count = Stack%Count + 1
       CALL SetErrorFlags(Stack,.FALSE.,0)
    ENDIF

  END SUBROUTINE add_to_stack

  ! (4) Remove the top node from the Stack (Pop)

  SUBROUTINE delete_from_stack(Stack,Node)
    IMPLICIT NONE
    TYPE(solution), OPTIONAL :: Node
    TYPE(StackT) :: Stack

    IF ( Stack%Last <= 0 ) THEN
       CALL SetErrorFlags(Stack,.TRUE.,-4)
       !IF ( PRESENT(Node) ) Node = HUGE(Node)
    ELSE
       IF ( PRESENT(Node) ) Node = select_from_stack(Stack)
       Stack%Last = Stack%Last - 1
       Stack%Count = Stack%Count - 1
       CALL SetErrorFlags(Stack,.FALSE.,0)
    ENDIF

  END SUBROUTINE delete_from_stack

  ! (5) Select the top node of the Stack for examination (Top/Peek)

  FUNCTION select_from_stack(Stack) RESULT(Node)
    IMPLICIT NONE
    TYPE(solution) :: Node
    TYPE(StackT) :: Stack

    IF ( Stack%Last >= 1 .AND. Stack%Last <= Stack%Size ) THEN
       Node = Stack%Key(Stack%Last)
       CALL SetErrorFlags(Stack,.FALSE.,0)
    ELSE
       !Node = HUGE(Node)
       CALL SetErrorFlags(Stack,.TRUE.,-5)
    ENDIF

  END FUNCTION select_from_stack

  ! (6) Convert Stack to an Array

  SUBROUTINE Stack2Array(Stack,Array)
    IMPLICIT NONE
    TYPE(solution), DIMENSION(1:) :: Array
    TYPE(StackT) :: Stack

    IF ( SIZE(Array) >= Stack%Last ) THEN
       array(1:Stack%Last) = Stack%Key(1:Stack%Last)
       CALL SetErrorFlags(Stack, .FALSE., 0)
    ELSE
       array = Stack%Key(1:SIZE(Array))
       CALL SetErrorFlags(Stack, .TRUE., -6)
    ENDIF

  END SUBROUTINE Stack2Array  

  ! (7) Check whether Stack has error flag on

  FUNCTION HasStackError(Stack) RESULT(Answer)
    IMPLICIT NONE
    LOGICAL :: Answer
    TYPE(StackT) :: Stack

    Answer = Stack%Error

  END FUNCTION HasStackError

  ! (8) Check the error number of Stack

  FUNCTION IsStackError(Stack) RESULT(n)
    IMPLICIT NONE
    INTEGER :: n
    TYPE(StackT) :: Stack
    
    n = Stack%ErrorNr

  END FUNCTION IsStackError


  ! (9) Check whether Stack is empty

  FUNCTION IsStackEmpty(Stack) RESULT(Answer)
    IMPLICIT NONE
    LOGICAL :: Answer
    TYPE(StackT) :: Stack
    
    Answer = ( Stack%Last == 0 )

  END FUNCTION IsStackEmpty


  ! (10) Check whether Stack is ful

  FUNCTION IsStackFul(Stack) RESULT(Answer)
    IMPLICIT NONE
    LOGICAL :: Answer
    TYPE(StackT) :: Stack

    Answer = ( Stack%Last == Stack%Size )

  END FUNCTION IsStackFul

  ! (11) How many Nodes populate the Stack?

  FUNCTION StackCardinality(Stack) RESULT(n)
    IMPLICIT NONE 
    INTEGER :: n
    TYPE(StackT) :: Stack

    n = Stack%Last

  END FUNCTION StackCardinality

  ! (12) Clear Stack and prepare it for new use

  SUBROUTINE ClearStack(Stack)
    IMPLICIT NONE
    TYPE(StackT) :: Stack

    Stack%Last = 0
    
  END SUBROUTINE ClearStack
  
    ! (13) get the number of element in this stack

  FUNCTION getCount(Stack) RESULT(stackSize)
    IMPLICIT NONE
    TYPE(StackT), intent(in) :: Stack
    integer :: stackSize
    
    stackSize = Stack%Count

  END FUNCTION getCount

  !------------ Sevice Routines -------------------!

  ! Set error flags and error numbers of Stack

  SUBROUTINE SetErrorFlags(Stack, Happened, Value)
    IMPLICIT NONE
    LOGICAL :: Happened
    INTEGER :: Value
    TYPE(StackT) :: Stack
  
    Stack%Error = Happened
    Stack%ErrorNr = Value

  END SUBROUTINE SetErrorFlags

END MODULE StackObject



    

    
  
  