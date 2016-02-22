from datetime import datetime
from fatman import db
from fatman.models import Task, Test, TestStructure, Method, TaskStatus

def main():
    """Given a list of test names and a method, add new tasks to the task list"""

    desired_tests = ["deltatest_H", "deltatest_Ne"]

    desired_method = Method.get(Method.id == 2)
    status_new = TaskStatus.get(TaskStatus.name == "new")

    
    all_teststructures = TestStructure.select()

    for x in all_teststructures:
        if x.test.name in desired_tests:
            Task.create(structure = x.structure, 
                        method    = desired_method, 
                        status    = status_new,
                        ctime     = datetime.now(),
                        mtime     = datetime.now(),
                        machine   = "-")


if __name__=="__main__":
    main()
