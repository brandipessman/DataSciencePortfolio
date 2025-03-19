SELECT posttest.recipient_id, posttest.question_id, question, answer_type, definition_1, definition_2, definition_3, definition_4, definition_5, answer
FROM training.posttest
JOIN training.questions
ON posttest.question_id = questions.question_id
WHERE posttest.question_id LIKE 'post26%';