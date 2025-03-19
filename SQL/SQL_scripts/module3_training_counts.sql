SELECT iteration, module, module3.question_id, definition_1, definition_2, definition_3, definition_4, definition_5, question, answer, COUNT(module3.recipient_id) AS response
FROM training.module3
JOIN training.participants
ON module3.recipient_id = participants.recipient_id
JOIN training.questions
ON module3.question_id = questions.question_id
WHERE module3.question_id NOT LIKE 'eval3_%_post_5%'
GROUP BY iteration, module, module3.question_id, definition_1, definition_2, definition_3, definition_4, definition_5, question, answer
ORDER BY iteration, question_id, answer;