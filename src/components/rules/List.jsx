import { RuleItem } from './Item'

export const RuleList = props => (
  <div>
    <div>
      <button type="button" onClick={props.handleAdd}>Add</button>
    </div>
    <ul>
      {props.rules.map((rule, index) => (
        <RuleItem
          key={rule.uuid}
          data={rule}
          domains={props.domains}
          rules={props.rules}
          handleRemove={props.handleRemove(index)}
          handleChange={props.handleChange(index)}
        />
      ))}
    </ul>
  </div>
)
