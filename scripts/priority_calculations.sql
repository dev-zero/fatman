CREATE VIEW resultwithouttestresult AS 
SELECT result.* FROM result JOIN task a  ON a.id = result.task_id 
    WHERE NOT EXISTS 
    (SELECT 1 FROM testresult b 
         WHERE a.method_id = b.method_id AND 
               a.test_id = b.test_id
    );
