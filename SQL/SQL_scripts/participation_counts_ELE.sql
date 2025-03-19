SELECT iteration, COUNT(answer) AS answer_count
FROM training.posttest
JOIN training.participants
ON posttest.recipient_id = participants.recipient_id
JOIN training.questions
ON posttest.question_id = questions.question_id
WHERE posttest.question_id ='participation5'
GROUP BY iteration;