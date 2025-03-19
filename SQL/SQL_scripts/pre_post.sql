select pretest.recipient_id, iteration, pretest.question_id, question, answer_type, definition_1, definition_2, definition_3, definition_4, definition_5, pretest.answer as pre_answer, posttest.answer as post_answer
FROM training.pretest
JOIN training.posttest
ON pretest.recipient_id = posttest.recipient_id AND pretest.question_id = posttest.question_id
JOIN training.questions
ON pretest.question_id = questions.question_id
JOIN training.participants
ON pretest.recipient_id = participants.recipient_id
WHERE pretest.answer IS NOT NULL AND posttest.answer IS NOT NULL
AND pretest.question_id NOT LIKE 'prepost3_%'
ORDER BY question_id;