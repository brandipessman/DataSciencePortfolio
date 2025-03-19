SELECT posttest.recipient_id, iteration, COUNT(answer) AS answer_count
FROM training.posttest
JOIN training.participants
ON posttest.recipient_id = participants.recipient_id
JOIN training.questions
ON posttest.question_id = questions.question_id
WHERE posttest.question_id LIKE 'participation%'
AND posttest.question_id NOT LIKE 'participation5'
GROUP BY posttest.recipient_id, iteration
ORDER BY posttest.recipient_id;