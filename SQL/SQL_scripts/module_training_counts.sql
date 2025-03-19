SELECT iteration, module, module1.question_id, definition_1, definition_2, definition_3, definition_4, definition_5, question, answer, COUNT(module1.recipient_id) AS response
FROM training.module1
JOIN training.participants
ON module1.recipient_id = participants.recipient_id
JOIN training.questions
ON module1.question_id = questions.question_id
GROUP BY iteration, module, module1.question_id, definition_1, definition_2, definition_3, definition_4, definition_5, question, answer
ORDER BY iteration, question_id, answer;