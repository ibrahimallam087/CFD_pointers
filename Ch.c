#include <stdio.h>
#include <stdlib.h>


typedef struct Node{
    int data;
    struct Node *next;

}Node;
void insertAtBeginning(Node **head,int value){
    Node * newNode = (Node*)malloc(sizeof(Node));
    if(newNode==NULL){
        printf("memory allocation failed!!\n");
        return;
    }

    newNode->data  = value;
    newNode->next = *head;
    *head = newNode;
}

void printList(Node *head){
    Node *temp = head;
    while(temp!= NULL){
        printf("%d->\n",temp->data);
        temp = temp->next;
    }
    printf("NULL\n");


}

int main(void){
    Node *head = NULL;

    insertAtBeginning(&head, 30);
    insertAtBeginning(&head, 20);
    insertAtBeginning(&head, 10);
    insertAtBeginning(&head, 5);

    printList(head); // Output: 5 → 10 → 20 → 30 → NULL

    return 0;





}
