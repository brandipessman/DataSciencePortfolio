SELECT posttest.recipient_id, iteration, posttest.question_id, question, definition_1, answer
FROM training.posttest
JOIN training.participants
ON posttest.recipient_id = participants.recipient_id
JOIN training.questions
ON posttest.question_id = questions.question_id
WHERE posttest.question_id LIKE 'participation%';